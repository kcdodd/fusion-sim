/*
    Copyright (C) 2015  Joulesmith Energy Technologies, LLC

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

define(['utilities', 'spindle'], function (util, spindle){
    "use strict";

    var exports = {};

    var N = function(num) {
        return num.toFixed(20);
    };

    var speed_of_light = 2.998e8; // m/s


    exports.makeCylindricalParticlePusher = function(spec) {
        util.validate_object(spec, {
            radius : 'number', // meters
            height : 'number',
            nr : 'number',
            nz : 'number',
            dt : 'number',
            nparticles : 'number',
            particle_mass : 'number',
            particle_charge : 'number',

        });

        // physical quantities
        var h = spec.particle_charge * spec.dt / (2 * spec.particle_mass);
        var factor_r = 1/spec.radius;
        var factor_z = 1/spec.height;

        var out = {};

        var i, j, k;

        // create canvas element for webgl to work on
        var canvas = document.createElement("CANVAS");
        canvas.id = "webgl_canvas_CylindricalParticlePusher";
        canvas.width = spec.nr;
        canvas.height = spec.nz;
        canvas.style.display = "none";

        document.body.appendChild(canvas);
        out.canvas = canvas;

        var webgl = util.webGL(canvas);

        webgl.enableFloatTexture();

        var vertex_positions = webgl.addVertexData([
            [-1, 1],
            [1, 1],
            [-1, -1],
            [-1, -1],
            [1, 1],
            [1, -1]
        ]);

        var render_vertex_positions = webgl.addVertexData([
            [-1, -1],
            [1, -1],
            [-1, 1],
            [-1, 1],
            [1, -1],
            [1, 1]
        ]);

        // texture coordinets for vertices
        var texture_coordinates = webgl.addVertexData([
            [0.0,  1.0],
            [1.0,  1.0],
            [0.0,  0.0],
            [0.0,  0.0],
            [1.0,  1.0],
            [1.0,  0.0]
        ]);

        var render_texture_coordinates = webgl.addVertexData([
            [0.0,  0.0],
            [1.0,  0.0],
            [0.0,  1.0],
            [0.0,  1.0],
            [1.0,  0.0],
            [1.0,  1.0]
        ]);

        //
        // The particles
        //

        var nparticles_h = spec.nparticles;
        var nparticles_w = spec.nparticles;
        var nparticles = nparticles_h * nparticles_w;

        var field_params = {
            width : spec.nr,
            height : spec.nz,
            useFloat : true
        };

        var particle_params = {
            width: nparticles_w,
            height: nparticles_h,
            useFloat: true
        };

        var position_arr = new Float32Array(4 * nparticles);

        var position_tex = webgl.addTextureArray({
            width: nparticles_w,
            height: nparticles_h,
            array: position_arr,
            useFloat: true
        });

        var velocity_arr = new Float32Array(4 * nparticles);

        var velocity_tex = webgl.addTextureArray({
            width: nparticles_w,
            height: nparticles_h,
            array: velocity_arr,
            useFloat: true
        });

        // source of entropy for pseudo-random numbers
        var n_entropy = 1024;
        var entropy_arr = new Float32Array(4 * n_entropy * n_entropy);

        var random_bytes = new Uint32Array(4);

        // initialize entropy values
        for(i = 0; i < n_entropy * n_entropy; i++) {

            window.crypto.getRandomValues(random_bytes);
            entropy_arr[4 * i] = random_bytes[0] / 0xFFFFFFFF;
            entropy_arr[4 * i + 1] = random_bytes[1] / 0xFFFFFFFF;
            entropy_arr[4 * i + 2] = random_bytes[2] / 0xFFFFFFFF;
            entropy_arr[4 * i + 3] = random_bytes[3] / 0xFFFFFFFF;
        }

        var entropy_tex = webgl.addTextureArray({
            width : n_entropy,
            height : n_entropy,
            array : entropy_arr,
            useFloat : true
        });

        // time dependant random numbers per particle
        var rand_arr = new Float32Array(4 * nparticles);

        // initialize random values
        for(i = 0; i < nparticles; i++) {
            rand_arr[4 * i] = Math.random();
            rand_arr[4 * i + 1] = Math.random();
            rand_arr[4 * i + 2] = Math.random();
            rand_arr[4 * i + 3] = Math.random();
        }

        var rand_tex = webgl.addTextureArray({
            width : nparticles_w,
            height : nparticles_h,
            array : rand_arr,
            useFloat : true
        });

        //
        // The electromagnetic fields
        //

        var E = webgl.addFrameBuffer(field_params);

        var E_arr = new Float32Array(4 * spec.nr * spec.nz);

        var E_tex = webgl.addTextureArray({
            width: spec.nr,
            height: spec.nz,
            array: E_arr,
            useFloat : true
        });

        var B = webgl.addFrameBuffer(field_params);

        var B_arr = new Float32Array(4 * spec.nr * spec.nz);

        var B_tex = webgl.addTextureArray({
            width: spec.nr,
            height: spec.nz,
            array: B_arr,
            useFloat : true
        });


        //
        // Mask for determining if particle is lost
        //

        var sink_mask = webgl.addFrameBuffer(field_params);

        var sink_mask_arr = new Float32Array(4 * spec.nr * spec.nz);

        var sink_mask_tex = webgl.addTextureArray({
            width: spec.nr,
            height: spec.nz,
            array: sink_mask_arr,
            useFloat: true
        });

        //
        // Inverse cdf in 2D for generating particle positions
        //

        var inv_cdf = webgl.addFrameBuffer({
            width: 512,
            height: 512,
            useFloat: true
        });

        var inv_cdf_arr = new Float32Array(4 * 512 * 512);

        var inv_cdf_tex = webgl.addTextureArray({
            width: 512,
            height: 512,
            array: inv_cdf_arr,
            useFloat : true
        });

        //
        // General shader to compute a texture frame buffer.
        //
        var precompute_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];

            return src_arr.join('\n');
        }; // precompute_vert()

        var avg_frag = function() {
            var src_arr = [
                "precision mediump float;",

                "uniform float u_ratio;",

                "uniform sampler2D u_next;",
                "uniform sampler2D u_avg;",
                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                    "vec4 next = texture2D(u_next, v_texCoord);",
                    "vec4 avg = texture2D(u_avg, v_texCoord);",

                    "gl_FragColor = u_ratio * next + (1.0 - u_ratio) * avg;",
                "}"
            ];

            return src_arr.join('\n');
        }; // avg_frag


        //
        // Program to calculate magnetic field from a current loop.
        //

        var B_loop_half = webgl.addFrameBuffer(field_params);

        var B_loop_tenth = webgl.addFrameBuffer(field_params);

        // ---------------------------------------------------------------------
        // computes general shape of field which can be scaled/translated for any loop.
        var programCurrentLoopShape = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_R;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        "float Bx = 0.0;",
                        "float Bz = 0.0;",
                        "float r = 0.0;",
                        "float factor = 0.0;",
                        "float cosine = 0.0;",
                        "float constant = u_R * 0.001 * 1.25663706e-6 / (4.0 * 3.14159265359);",

                        "for(float i = 0.0; i < 1000.0; i++){",

                            "cosine = cos(3.14159265359 * (i + 0.5) / 1000.0);",
                            "r = sqrt(u_R * u_R + v_texCoord.x * v_texCoord.x + v_texCoord.y * v_texCoord.y - 2.0 * v_texCoord.x * u_R * cosine);",
                            "factor = (r > 0.0) ? constant * (1.0) / (r * r * r) : 0.0;",

                            "Bx += v_texCoord.y * factor * cosine;",
                            "Bz += factor * (u_R - v_texCoord.x * cosine);",

                        "}",

                        "gl_FragColor = vec4(Bx, 0.0, Bz, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_R" : 0.5 // compute one with loop in middle for positions close to axis
        }).draw({
            triangles : 6,
            target : B_loop_half
        }).set({
            "u_R" : 0.1 // compute one with small radius to be scaled for large distances
        }).draw({
            triangles : 6,
            target : B_loop_tenth
        });

        // ---------------------------------------------------------------------
        // program for computing field from current loop
        var programCurrentLoop = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_R;",
                    "uniform float u_Z;",
                    "uniform float u_I;",

                    "uniform sampler2D u_shape_half;",
                    "uniform sampler2D u_shape_tenth;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        "float a = v_texCoord.x / u_R;",
                        "float b = (v_texCoord.y - u_Z) / u_R;",

                        "vec4 field;",
                        "if(a > 2.0 || b > 2.0){",
                            "field = u_I * vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "field = u_I * vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",

                        "gl_FragColor = field;",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_shape_half" : B_loop_half,
            "u_shape_tenth" : B_loop_tenth
        });

        // ---------------------------------------------------------------------
        var programCurrentZ = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_I;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, u_I * 1.25663706e-6 / (2.0 * 3.14159265359 * v_texCoord.x), 0.0, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        var programBZ = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_Bz;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, 0.0, u_Bz, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        var programBTheta = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_Btheta;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor += vec4(0.0, u_Btheta, 0.0, 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        var programBMag = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_B;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float mag = length(B);",
                        "vec3 dir = B / mag;",
                        "gl_FragColor =  vec4(mag * abs (min(0.0, dir.z)), mag * dir.x, mag * abs (max(0.0, dir.z)), 1.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_B" : B
        });

        //
        // programs for pre-calculations of the fields
        //

        var R1 = webgl.addFrameBuffer(field_params);
        var R2 = webgl.addFrameBuffer(field_params);
        var R3 = webgl.addFrameBuffer(field_params);
        var A = webgl.addFrameBuffer(field_params);


        // ---------------------------------------------------------------------
        var programPre1 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R11 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.x * B.x;",
                        "float R12 = factor * u_h * (B.z + u_h * B.x * B.y);",
                        "float R13 = factor * u_h * (-B.y + u_h * B.x * B.z);",
                        "vec3 R1 = vec3(R11, R12, R13 * " + N(factor_r/factor_z) + ");",
                        //"vec3 R1 = vec3(R11 , R12, R13);",

                        // R1
                        "gl_FragColor = vec4(R1, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_B" : B,
            "u_h" : h
        });

        // ---------------------------------------------------------------------
        var programPre2 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R21 = factor * u_h * (-B.z + u_h * B.y * B.x);",
                        "float R22 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.y * B.y;",
                        "float R23 = factor * u_h * (B.x + u_h * B.y * B.z);",
                        "vec3 R2 = vec3(R21, R22, R23 * " + N(factor_r/factor_z) + ");",
                        //"vec3 R2 = vec3(R21 , R22, R23);",

                        // R2
                        "gl_FragColor = vec4(R2, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_B" : B,
            "u_h" : h
        });


        // ---------------------------------------------------------------------
        var programPre3 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_B;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0 / (1.0 + hB2);",

                        "float R31 = factor * u_h * (B.y + u_h * B.z * B.x);",
                        "float R32 = factor * u_h * (-B.x + u_h * B.z * B.y);",
                        "float R33 = (1.0 - hB2 * factor) + factor * u_h * u_h * B.z * B.z;",
                        "vec3 R3 = vec3(R31 * " + N(factor_z/factor_r) + ", R32 * " + N(factor_z/factor_r) + ", R33);",
                        //"vec3 R3 = vec3(R31 , R32, R33);",

                        // R2
                        "gl_FragColor = vec4(R3, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_B" : B,
            "u_h" : h
        });


        // ---------------------------------------------------------------------
        var programPreA = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_h;",

                    "uniform sampler2D u_E;",
                    "uniform sampler2D u_B;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec3 B = texture2D(u_B, v_texCoord).xyz;",
                        "vec3 E = texture2D(u_E, v_texCoord).xyz;",
                        "float Bmag = length(B);",
                        "float hB2 = u_h * u_h * Bmag * Bmag;",
                        "float factor = 2.0/(1.0 + hB2);",
                        "vec3 A = (u_h * (2.0 - hB2 * factor) * E + u_h * u_h * factor * (cross(E,B) + u_h * dot(E,B)))/2.998e8;",
                        // divide by speed of light because of normalized velocity
                        "gl_FragColor = vec4(A.x * " + N(factor_r) + ", A.y * " + N(factor_r) + ", A.z * " + N(factor_z) + ", 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_B" : B,
            "u_E" : E,
            "u_h" : h
        });


        //
        // step programs in leap-frog
        //

        var position_A = webgl.addFrameBuffer(particle_params);
        var velocity_A = webgl.addFrameBuffer(particle_params);
        var position_B = webgl.addFrameBuffer(particle_params);
        var velocity_B = webgl.addFrameBuffer(particle_params);

        var rand_A = webgl.addFrameBuffer(particle_params);
        var rand_B = webgl.addFrameBuffer(particle_params);

        // ---------------------------------------------------------------------
        var step_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                "varying vec2 v_texCoord;",

                "void main() {",
                "    gl_Position = vec4(a_position, 0, 1);",

                "    v_texCoord = a_texCoord;",
                "}",
            ];

            return src_arr.join('\n');
        } // step_vert()

        // ---------------------------------------------------------------------
        var step_position_frag = function() {
            var src_arr = [
                "precision highp float;",

                "uniform float u_step_factor;",

                // current phase
                "uniform sampler2D u_position;",
                "uniform sampler2D u_velocity;",

                // current random
                "uniform sampler2D u_rand;",

                // sources and sinks
                "uniform sampler2D u_sink;",
                "uniform sampler2D u_inv_cdf;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",

                    "vec4 next_position = texture2D(u_position, v_texCoord) + u_step_factor * texture2D(u_velocity, v_texCoord);",
                    "float r = sqrt(next_position.x * next_position.x + next_position.y * next_position.y);",
                    "vec4 rand = texture2D(u_rand, v_texCoord);",
                    "vec2 new_position = texture2D(u_inv_cdf, vec2(rand.x, rand.y)).xy;",

                    "gl_FragColor = texture2D(u_sink, vec2(r, next_position.z)).r > 0.5 ? vec4(next_position.xyz, 1.0) : vec4(new_position.x, 0.0, new_position.y, 0.0);",

                "}"
            ];


            return src_arr.join('\n');
        }; // step_position_frag()

        // ---------------------------------------------------------------------
        var step_velocity_frag = function() {
            var src_arr = [
                "precision highp float;",

                // current phase
                "uniform sampler2D u_position;",
                "uniform sampler2D u_velocity;",

                // current random
                "uniform sampler2D u_rand;",

                // fields
                "uniform sampler2D u_R_1;",
                "uniform sampler2D u_R_2;",
                "uniform sampler2D u_R_3;",
                "uniform sampler2D u_A;",

                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                    "vec4 rand = texture2D(u_rand, v_texCoord);",
                    "vec4 position4 = texture2D(u_position, v_texCoord);",
                    "vec3 position = position4.xyz;",
                    "vec3 velocity = texture2D(u_velocity, v_texCoord).xyz;",

                    "float r = sqrt(position.x * position.x + position.y * position.y);",
                    "vec2 direction = vec2(position.x/r, position.y/r);",

                    "float vr = velocity.x * direction.x + velocity.y * direction.y;",
                    "float va = velocity.y * direction.x - velocity.x * direction.y;",
                    "vec2 cylindrical_position = vec2(r, position.z);",
                    "vec3 cylindrical_velocity = vec3(vr, va, velocity.z);",

                    "vec3 R1 = texture2D(u_R_1, cylindrical_position).xyz;",
                    "vec3 R2 = texture2D(u_R_2, cylindrical_position).xyz;",
                    "vec3 R3 = texture2D(u_R_3, cylindrical_position).xyz;",
                    "vec3 A = texture2D(u_A, cylindrical_position).xyz;",

                    "cylindrical_velocity = vec3(dot(R1, cylindrical_velocity), dot(R2, cylindrical_velocity), dot(R3, cylindrical_velocity)) + A;",
                    "vec3 next_velocity = vec3(cylindrical_velocity.x * direction.x - cylindrical_velocity.y * direction.y, cylindrical_velocity.x * direction.y + cylindrical_velocity.y * direction.x, cylindrical_velocity.z);",

                    // position4.a === 0 indicates that the particle was just generated and needs a new velocity
                    "gl_FragColor = position4.a > 0.5 ? vec4(next_velocity, 1.0) : 0.001 * vec4(2.0 * rand.x - 1.0, 2.0 * rand.y - 1.0, 2.0 * rand.z - 1.0, 1.0);",
                "}"
            ];


            return src_arr.join('\n');
        }; // step_velocity_frag()


        // ---------------------------------------------------------------------
        // computes value of rand_B
        var programStepRandB = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    // source of entropy
                    "uniform sampler2D u_entropy;",

                    // current random
                    "uniform sampler2D u_rand;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        "vec4 rand = texture2D(u_rand, v_texCoord);",
                        "vec2 x = rand.zw;",
                        "vec4 s = texture2D(u_entropy, x);",

                        "x = 0.999 * x + 0.001 * s.zw;",
                        "vec2 m = rand.xy + s.xy;",
                        //
                        "gl_FragColor = vec4((m.x > 1.0) ? m.x - 1.0 : m.x, (m.y > 1.0) ? m.y - 1.0 : m.y, 4.0 * x * (1.0 - x));",


                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_entropy" : entropy_tex,
            "u_rand" : rand_A
        });

        // ---------------------------------------------------------------------
        // computes value of velocity_B
        var programStepVelocityB = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_velocity_frag()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_position" : position_A,
            "u_velocity" : velocity_A,
            "u_rand" : rand_A,
            "u_R_1" : R1,
            "u_R_2" : R2,
            "u_R_3" : R3,
            "u_A" : A
        });

        // ---------------------------------------------------------------------
        // computes value of position_B
        var programStepPositionB = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_position_frag()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_position" : position_A,
            "u_velocity" : velocity_B,
            "u_rand" : rand_A,
            "u_sink" : sink_mask,
            "u_inv_cdf" : inv_cdf,
            "u_step_factor" : spec.dt * speed_of_light
        });


        // ---------------------------------------------------------------------
        // computes value of rand_A
        var programStepRandA = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    // source of entropy
                    "uniform sampler2D u_entropy;",

                    // current random
                    "uniform sampler2D u_rand;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        "vec4 rand = texture2D(u_rand, v_texCoord);",
                        "vec2 x = rand.zw;",
                        "vec4 s = texture2D(u_entropy, x);",

                        "x = 0.999 * x + 0.001 * s.zw;",
                        "vec2 m = rand.xy + s.xy;",
                        //
                        "gl_FragColor = vec4((m.x > 1.0) ? m.x - 1.0 : m.x, (m.y > 1.0) ? m.y - 1.0 : m.y, 4.0 * x * (1.0 - x));",


                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_entropy" : entropy_tex,
            "u_rand" : rand_B
        });

        // ---------------------------------------------------------------------
        // computes value of velocity_A
        var programStepVelocityA = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_velocity_frag()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_position" : position_B,
            "u_velocity" : velocity_B,
            "u_rand" : rand_B,
            "u_R_1" : R1,
            "u_R_2" : R2,
            "u_R_3" : R3,
            "u_A" : A
        });

        // ---------------------------------------------------------------------
        // computes value of position_A
        var programStepPositionA = webgl.linkProgram({
            vertexShaderSource : step_vert(),
            fragmentShaderSource : step_position_frag()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_position" : position_B,
            "u_velocity" : velocity_A,
            "u_rand" : rand_B,
            "u_sink" : sink_mask,
            "u_inv_cdf" : inv_cdf,
            "u_step_factor" : spec.dt * speed_of_light
        });

        // ---------------------------------------------------------------------
        // program for rendering particles into a density and velocity function

        var moments01 = webgl.addFrameBuffer(field_params);


        // texture coordinates of particles
        var particles_texcoord_arr = [];

        for(j = 0; j < nparticles_h; j++) {
            for(i = 0; i < nparticles_w; i++) {
                particles_texcoord_arr[i + nparticles_w * j] = [(i+0.5)/(nparticles_w), (j+0.5)/(nparticles_h)];
            }
        }


        var particles_texcoord = webgl.addVertexData(particles_texcoord_arr);


        var nshape = 11;
        var shape_arr = new Float32Array(4 * nshape * nshape);

        var mid = (nshape-1)/2;

        var sum = 0;

        for(j = 0; j < nshape; j++) {
            for(i = 0; i < nshape; i++) {
                var d = Math.sqrt(Math.pow(i - mid, 2) + Math.pow(j - mid, 2));
                shape_arr[4*(i + nshape*j)] = Math.pow(Math.max(0.0, Math.cos(0.5*Math.PI*d/mid)), 2);
                sum += shape_arr[4*(i + nshape*j)];
            }
        }

        for(j = 0; j < nshape; j++) {
            for(i = 0; i < nshape; i++) {
                shape_arr[4*(i + nshape*j)] = shape_arr[4*(i + nshape*j)] / sum;
                shape_arr[4*(i + nshape*j)+1] = shape_arr[4*(i + nshape*j)];
                shape_arr[4*(i + nshape*j)+2] = shape_arr[4*(i + nshape*j)];
                shape_arr[4*(i + nshape*j)+3] = shape_arr[4*(i + nshape*j)];
            }
        }

        var shape_tex = webgl.addTextureArray({
            width: nshape,
            height: nshape,
            array: shape_arr,
            useFloat: true
        });

        var programMoments01 = webgl.linkProgram({
            vertexShaderSource : (function() {
                var src_arr = [
                    "attribute vec2 a_particleTexCoord;",

                    "uniform sampler2D u_position;",
                    "uniform sampler2D u_velocity;",

                    "uniform float u_pointsize;",

                    "uniform float u_weight;",
                    "varying vec4 v_color;",

                    "void main() {",
                        "vec3 position = texture2D(u_position, a_particleTexCoord).xyz;",
                        "vec3 velocity = texture2D(u_velocity, a_particleTexCoord).xyz;",
                        "float r = sqrt(position.x * position.x + position.y * position.y);",
                        "gl_Position = vec4(2.0*r-1.0, 2.0*position.z-1.0, 0, 1.0);",
                        //"gl_Position = vec4(position.x, position.y, 0, 1.0);",
                        "gl_PointSize = u_pointsize;",

                        "vec2 direction = vec2(position.x/r, position.y/r);",

                        "float v = length(velocity);",
                        "float vr = velocity.x * direction.x + velocity.y * direction.y;",
                        "float va = velocity.y * direction.x - velocity.x * direction.y;",
                        "v_color = 0.001 * vec4(vr, va, velocity.z, 1.0);",
                        //"v_color = vec4(1.0, 1.0, 1.0, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_shape;",

                    "varying vec4 v_color;",

                    "void main() {",
                        "gl_FragColor = v_color * texture2D(u_shape, gl_PointCoord);",
                    "}"
                ];


                return src_arr.join('\n');
            })()
        }).set({
            "a_particleTexCoord" : particles_texcoord,
            "u_position" : position_A,
            "u_velocity" : velocity_A,
            "u_shape" : shape_tex,
            "u_pointsize" : nshape
        });

        // ---------------------------------------------------------------------
        // program for normalizing average of moment
        //
        var moments01_norm = webgl.addFrameBuffer(field_params);

        var programNormalizeMoments01 = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_moments01;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec4 M01 = texture2D(u_moments01, v_texCoord);",
                        //"gl_FragColor = vec4(M01.xyz/M01.w, M01.w);",
                        "M01 = (M01.a > 0.0) ? vec4(M01.r / M01.a, M01.g / M01.a, M01.b / M01.a, M01.a) : vec4(0.0, 0.0, 0.0, 0.0);",
                        "gl_FragColor = 1000.0 * M01 * 0.5 / v_texCoord.x;",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_moments01" : moments01
        });

        // ---------------------------------------------------------------------
        // program for normalizing average of moment
        //
        var moments01_avgA = webgl.addFrameBuffer(field_params);
        var moments01_avgB = webgl.addFrameBuffer(field_params);

        // ---------------------------------------------------------------------
        var programAvgMoments = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : avg_frag()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_next" : moments01_norm,
            "u_avg" : moments01_avgB,
            "u_ratio" : 0.01
        });


        // ---------------------------------------------------------------------
        // program for rendering density as a color
        //
        var programDensity = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_moments01;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec4 M01 = texture2D(u_moments01, v_texCoord);",
                        "float v = length(M01.xyz);",
                        "float va = v > 0.0 ? M01.y / v : 0.0;",
                        //"gl_FragColor = vec4(M01.a * (1.0 + min(0.0, va)), M01.a * (1.0 - abs(va)), M01.a * (1.0 - max(0.0, va)), 1.0);",
                        "gl_FragColor = 0.5 * vec4(M01.a, M01.a, M01.a, 1.0);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_moments01" : moments01_avgA,
            //"u_moments01" : moments01_norm
        });

        // ---------------------------------------------------------------------
        // program for setting value of particles
        //
        var programSet = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_value;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "gl_FragColor = texture2D(u_value, v_texCoord);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        // initialize random values

        programSet.set({
            "u_value" : rand_tex
        }).draw({
            triangles : 6,
            target : rand_A
        });

        //
        // Output interface for simulation programs
        //

        out.set = function(value) {

            if (value.E) {
                for(i = 0; i < spec.nr; i++) {
                    for(j = 0; j < spec.nz; j++){
                        E_arr[4*(i + j*spec.nr)] = value.E[i][j][0];
                        E_arr[4*(i + j*spec.nr)+1] = value.E[i][j][1];
                        E_arr[4*(i + j*spec.nr)+2] = value.E[i][j][2];
                        E_arr[4*(i + j*spec.nr)+3] = 1.0;
                    }
                }

                E_tex.update();

                programSet.set({"u_value" : E_tex});

                programSet.draw({
                    triangles : 6,
                    target : E
                });
            }

            if (value.B) {
                for(i = 0; i < spec.nr; i++) {
                    for(j = 0; j < spec.nz; j++){
                        B_arr[4*(i + j*spec.nr)] = value.B[i][j][0];
                        B_arr[4*(i + j*spec.nr)+1] = value.B[i][j][1];
                        B_arr[4*(i + j*spec.nr)+2] = value.B[i][j][2];
                        B_arr[4*(i + j*spec.nr)+3] = 1.0;
                    }
                }

                B_tex.update();

                programSet.set({"u_value" : B_tex});

                programSet.draw({
                    triangles : 6,
                    target : B
                });
            }

            if (value.position) {

                for(i = 0; i < nparticles; i++) {
                    position_arr[4*i] = value.position[i][0] * factor_r;
                    position_arr[4*i+1] = value.position[i][1] * factor_r;
                    position_arr[4*i+2] = value.position[i][2] * factor_z;
                    position_arr[4*i+3] = 1.0;
                }

                // send values to card
                position_tex.update();
                programSet.set({"u_value" : position_tex});

                programSet.draw({
                    triangles : 6,
                    target : position_A
                });

                programSet.draw({
                    triangles : 6,
                    target : position_B
                });
            }

            if (value.velocity) {

                for(i = 0; i < nparticles; i++) {
                    velocity_arr[4*i] = value.velocity[i][0] * factor_r;
                    velocity_arr[4*i+1] = value.velocity[i][1] * factor_r;
                    velocity_arr[4*i+2] = value.velocity[i][2] * factor_z;
                    velocity_arr[4*i+3] = 1.0;
                }

                // send values to card
                velocity_tex.update();
                programSet.set({"u_value" : velocity_tex});

                programSet.draw({
                    triangles : 6,
                    target : velocity_A
                });
                programSet.draw({
                    triangles : 6,
                    target : velocity_B
                });
            }

            if (value.sink_mask) {
                for(i = 0; i < spec.nr; i++) {
                    for(j = 0; j < spec.nz; j++){
                        sink_mask_arr[4*(i + j*spec.nr)] = value.sink_mask[i][j];
                    }
                }

                sink_mask_tex.update();
                programSet.set({"u_value" : sink_mask_tex});

                programSet.draw({
                    triangles : 6,
                    target : sink_mask
                });
            }


            if (value.source_pdf) {

                // a pdf cannot be used directly to generation particle positions


                var cdf_y = [];
                var cdf_x = [];
                var sum_x = 0;

                for(i = 0; i < value.source_pdf.length; i++) {
                    cdf_y[i] = [];
                    var sum_y = 0;

                    for(j = 0; j < value.source_pdf[0].length; j++) {
                        sum_y += value.source_pdf[i][j];
                        cdf_y[i][j] = sum_y;
                    }

                    for(j = 0; j < value.source_pdf[0].length; j++) {
                        cdf_y[i][j] /= sum_y;
                    }

                    sum_x += sum_y;
                    cdf_x[i] = sum_x;
                }

                for(i = 0; i < value.source_pdf.length; i++) {
                    cdf_x[i] /= sum_x;
                }

                var inverse_cdf_x = function(f) {
                    if (f < 0 || f > 1) {
                        throw new Error("function out of range");
                    }

                    var i = 0;

                    while(cdf_x[i] < f) {
                        i++;
                    }

                    if (i === 0) {
                        return (f/cdf_x[0])/value.source_pdf.length;
                    }

                    return (i + (f - cdf_x[i-1])/(cdf_x[i] - cdf_x[i-1]))/value.source_pdf.length;
                };

                var inverse_cdf_y = function(x, f) {
                    var i = Math.min(value.source_pdf.length-1, Math.floor(x * value.source_pdf.length));

                    var j = 0;

                    while(cdf_y[i][j] < f) {
                        j++;

                    }

                    if (j === 0) {
                        return (f/cdf_y[i][0])/value.source_pdf[0].length;
                    }

                    return (j + (f - cdf_y[i][j-1])/(cdf_y[i][j] - cdf_y[i][j-1]))/value.source_pdf[0].length;
                };

                for(i = 0; i < 512; i++) {
                    var f1 = i/511;

                    for(j = 0; j < 512; j++) {
                        var f2 = j/511;
                        var x = inverse_cdf_x(f1);
                        var y = inverse_cdf_y(x, f2);

                        inv_cdf_arr[4 * (i + j * 512)] = x;
                        inv_cdf_arr[4 * (i + j * 512) + 1] = y;
                    }
                }

                inv_cdf_tex.update();
                programSet.set({"u_value" : inv_cdf_tex});

                programSet.draw({
                    triangles : 6,
                    target : inv_cdf
                });

            }
        };

        out.addCurrentLoop = function(r, z, I) {

            programCurrentLoop.set({
                "u_R" : r * factor_r,
                "u_Z" : z * factor_z,
                "u_I" : I
            }).draw({
                triangles : 6,
                target : B,
                blend : ['ONE', 'ONE']
            });
        };

        /**
            Solves the boundary conditions for a perfect conductor in center of
            a spindle cusp magnetic field.
        */
        out.addSpindleCuspPlasmaField = function(r, B_c, beta_c) {

            var currents = spindle.makeSpindleCuspPlasmaField({
                radius : spec.radius, // meters
                height : spec.height,
                nr : spec.nr,
                nz : spec.nz,
                webgl : webgl
            });
        };

        out.addCurrentZ = function(I) {

            programCurrentZ.set({
                "u_I" : I
            }).draw({
                triangles : 6,
                target : B,
                blend : ['ONE', 'ONE']
            });
        };

        out.addBZ = function(Bz) {

            programBZ.set({
                "u_Bz" : Bz
            }).draw({
                triangles : 6,
                target : B,
                blend : ['ONE', 'ONE']
            });
        };

        out.addBTheta = function(Btheta) {

            programBTheta.set({
                "u_Btheta" : Btheta
            }).draw({
                triangles : 6,
                target : B,
                blend : ['ONE', 'ONE']
            });
        };

        out.precalc = function () {

            programPre1.draw({
                triangles : 6,
                target : R1
            });

            programPre2.draw({
                triangles : 6,
                target : R2
            });

            programPre3.draw({
                triangles : 6,
                target : R3
            });

            programPreA.draw({
                triangles : 6,
                target : A
            });
        };

        out.step = function () {

            programStepRandB.draw({
                triangles : 6,
                target : rand_B
            });

            programStepVelocityB.draw({
                triangles : 6,
                target : velocity_B
            });

            programStepPositionB.draw({
                triangles : 6,
                target : position_B
            });


            programStepRandA.draw({
                triangles : 6,
                target : rand_A
            });

            programStepVelocityA.draw({
                triangles : 6,
                target : velocity_A
            });

            programStepPositionA.draw({
                triangles : 6,
                target : position_A
            });

        };

        out.density = function() {

            programMoments01.draw({
                points : nparticles,
                target : moments01,
                clear_color : [0,0,0,0],
                blend : ['ONE', 'ONE']
            });

            programNormalizeMoments01.draw({
                triangles : 6,
                target : moments01_norm
            });

            programAvgMoments.draw({
                triangles : 6,
                target : moments01_avgA
            });

            programSet.set({
                "u_value" : moments01_avgA
            }).draw({
                triangles : 6,
                target : moments01_avgB
            });

            programBMag.draw({
                triangles : 6,
            });

            programDensity.draw({
                triangles : 6,
                blend : ['SRC_ALPHA', 'ONE']
            });


/*
            programStepRandB.draw({
                triangles : 6,
                target : rand_B,
                //clear_color : [0,0,0,0],
            });

            programStepRandA.draw({
                triangles : 6,
                target : rand_A,
                //clear_color : [0,0,0,0],
            });
*/
/*
            programSet.set({"u_value" : rand_A});
            programSet.draw({
                triangles : 6
            });
*/
        };

        return out;
    };

    return exports;
});
