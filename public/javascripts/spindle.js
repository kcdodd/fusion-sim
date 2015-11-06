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

define(['utilities', 'matrix_webgl'], function (util, matrix_webgl){
    "use strict";

    var exports = {};

    var N = function(num) {
        return num.toFixed(20);
    };

    /**
        Solves the boundary conditions for a perfect conductor in center of
        a spindle cusp magnetic field.
    */
    exports.makeSpindleCuspPlasmaField = function(spec) {
        util.validate_object(spec, {
            radius : 'number', // meters
            height : 'number',
            nr : 'number',
            nz : 'number',

            webgl : [,'object']
        });

        var factor_r = 1/spec.radius;
        var factor_z = 1/spec.height;

        var webgl;

        if (spec.webgl) {
            webgl = spec.webgl;
        }else{
            // create canvas element for webgl to work on
            var canvas = document.createElement("CANVAS");
            canvas.id = "webgl_canvas_SpindleCuspPlasmaField";
            canvas.width = spec.nr;
            canvas.height = spec.nz;
            canvas.style.display = "none";

            document.body.appendChild(canvas);
            out.canvas = canvas;

            var webgl = util.webGL(canvas);
        }

        webgl.enableFloatTexture();

        var eq = matrix_webgl.makeSORIterative({
            n_power : 3,
            webgl : webgl
        });

        var field_params = {
            width : spec.nr,
            height : spec.nz,
            useFloat : true
        };


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


        var n_loops = eq.vec_length;



        // these are n^2 redundant because this is used to make the matrix
        // and 4x because the red channel is only one used.
        // current loop positions
        var loops_arr_A = new Float32Array(4 * n_loops * n_loops);

        // positions to evaluate field
        var points_arr_A = new Float32Array(4 * n_loops * n_loops);

        // normals at positions to evaluate field
        var normals_arr_A = new Float32Array(4 * n_loops * n_loops);

        // same as above, but for computing the vector instead of matrix
        var pointsX_arr_b = new Float32Array(n_loops);
        var pointsZ_arr_b = new Float32Array(n_loops);

        var normalsX_arr_b = new Float32Array(n_loops);
        var normalsZ_arr_b = new Float32Array(n_loops);

        // this is the initial guess for the current distribution
        var currents_arr = new Float32Array(n_loops);


        var a = 0.4;
        var R = spec.radius * Math.sqrt(1 + a * a);
        // side angle to start arc
        var alpha =  Math.atan(a);

        // total starting angle
        var theta = alpha + Math.PI;
        // total angle of arc
        var arc = 0.5 * Math.PI - 2.0 * alpha;

        for(var point = 0; point < n_loops; point++) {

            var phi = (point + 0.5) * arc / 1000 + theta;

            var x = R * Math.cos(-phi) + spec.radius; // point.x
            var z = R * Math.sin(-phi); // point.z
            var n_x = - Math.cos(-phi); // normal.x
            var n_z = - Math.sin(-phi); // normal.z

            currents_arr[point] = 1.0;

            pointsX_arr_b[point] = x; // point.x
            pointsZ_arr_b[point] = z; // point.z

            normalsX_arr_b[point] = n_x; // normal.x
            normalsZ_arr_b[point] = n_z; // normal.z

            for(var loop = 0; loop < n_loops; loop++) {
                var phi_0 = (loop) * arc / 1000 + theta;
                var phi_1 = (loop + 1) * arc / 1000 + theta;

                var x0 = R * Math.cos(-phi_0) + spec.radius; // loop0.x
                var z0 = R * Math.sin(-phi_0); // loop0.z

                var x1 = R * Math.cos(-phi_1) + spec.radius; // loop1.x
                var z1 = R * Math.sin(-phi_1); // loop1.z

                var offset = 4 * (loop + point * n_loops);

                points_arr_A[offset] = x; // point.x
                points_arr_A[offset + 2] = z; // point.z

                normals_arr_A[offset] = n_x; // normal.x
                normals_arr_A[offset + 2] = n_z; // normal.z

                loops_arr_A[offset] = x0; // loop0.x
                loops_arr_A[offset + 1] = z0; // loop0.z
                loops_arr_A[offset + 2] = x1; // loop1.x
                loops_arr_A[offset + 3] = z1; // loop1.z
            }
        }

        // add the arrays to webgl to use as textures
        var loops_tex_A = webgl.addTextureArray({
            width: n_loops,
            height: n_loops,
            array: loops_arr_A,
            useFloat : true
        });

        var points_tex_A = webgl.addTextureArray({
            width: n_loops,
            height: n_loops,
            array: points_arr_A,
            useFloat : true
        });

        var normals_tex_A = webgl.addTextureArray({
            width: n_loops,
            height: n_loops,
            array: normals_arr_A,
            useFloat : true
        });

        var pointsX_tex_b = webgl.addTextureArray({
            width: eq.vec_height,
            height: eq.vec_height,
            array: pointsX_arr_b,
            useFloat : true
        });

        var pointsZ_tex_b = webgl.addTextureArray({
            width: eq.vec_height,
            height: eq.vec_height,
            array: pointsZ_arr_b,
            useFloat : true
        });

        var normalsX_tex_b = webgl.addTextureArray({
            width: eq.vec_height,
            height: eq.vec_height,
            array: normalsX_arr_b,
            useFloat : true
        });

        var normalsZ_tex_b = webgl.addTextureArray({
            width: eq.vec_height,
            height: eq.vec_height,
            array: normalsZ_arr_b,
            useFloat : true
        });

        //
        // General shader to compute a texture frame buffer.
        //
        var render_vert = function() {
            var src_arr = [
                "attribute vec2 a_position;",
                "attribute vec2 a_texCoord;",
                "varying vec2 v_texCoord;",

                "void main() {",
                    "gl_Position = vec4(a_position, 0, 1);",

                    "v_texCoord = a_texCoord;",
                "}",
            ];

            return src_arr.join('\n');
        }; // render_vert()


        //
        // Program to calculate magnetic field from a current loop.
        //

        var B_loop_half = webgl.addFrameBuffer(field_params);

        var B_loop_tenth = webgl.addFrameBuffer(field_params);

        // ---------------------------------------------------------------------
        // computes general shape of field which can be scaled/translated for any loop.
        var programCurrentLoopShape = webgl.linkProgram({
            vertexShaderSource : render_vert(),
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
        var programLoopField = webgl.linkProgram({
            vertexShaderSource : precompute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_I",

                    "uniform float u_R0;",
                    "uniform float u_Z0;",

                    "uniform float u_R1;",
                    "uniform float u_Z1;",

                    "uniform sampler2D u_shape_half;",
                    "uniform sampler2D u_shape_tenth;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        // loop 0
                        "float a = v_texCoord.x / u_R0;",
                        "float b = (v_texCoord.y - u_Z0) / u_R0;",

                        "vec4 field = vec4(0.0);",

                        "vec4 tmp = vec4(0.0);",
                        "if(a > 2.0 || b > 2.0){",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",
                        "field += tmp;",

                        // loop 0 reflected vertically
                        "a = v_texCoord.x / u_R0;",
                        "b = (v_texCoord.y - (1.0 - u_Z0)) / u_R0;",

                        "if(a > 2.0 || b > 2.0){",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",
                        "field += tmp;",

                        // loop 1
                        "a = v_texCoord.x / u_R1;",
                        "b = (v_texCoord.y - u_Z1) / u_R1;",

                        "if(a > 2.0 || b > 2.0){",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",
                        "field += tmp;",

                        // loop 1 reflected vertically
                        "a = v_texCoord.x / u_R1;",
                        "b = (v_texCoord.y - (1.0 - u_Z1)) / u_R1;",

                        "if(a > 2.0 || b > 2.0){",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",
                        "field += tmp;",

                        "gl_FragColor = u_I * field;",

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

        // compute vector b
        var boundary_b = webgl.addFrameBuffer({
            width: eq.vec_height,
            height: eq.vec_height,
            useFloat: true
        });

        // ---------------------------------------------------------------------
        // program for computing field normal contributions at points from fixed loops
        var programFixedLoop = webgl.linkProgram({
            vertexShaderSource : render_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform float u_R;",
                    "uniform float u_Z;",
                    "uniform float u_I;",

                    "uniform sampler2D u_shape_half;",
                    "uniform sampler2D u_shape_tenth;",

                    "uniform sampler2D u_pointsX;",
                    "uniform sampler2D u_pointsZ;",
                    "uniform sampler2D u_normalsX;",
                    "uniform sampler2D u_normalsZ;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec4 pointX = texture2D(u_pointsX, v_texCoord);",
                        "vec4 pointZ = texture2D(u_pointsZ, v_texCoord);",
                        "vec4 normalX = texture2D(u_normalsX, v_texCoord);",
                        "vec4 normalZ = texture2D(u_normalsZ, v_texCoord);",


                        "float a = pointX.r / u_R;",
                        "float b = (pointZ.r - u_Z) / u_R;",

                        "vec4 normal = u_I * vec4(sign(b) * normalX.r, 0.0, normalZ.r, 0.0);",
                        "float red;",

                        "if(a > 2.0 || b > 2.0){",
                            "red = dot(normal, texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0)));",
                        "}else{",
                            "red = dot(normal,texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0)));",
                        "}",

                        "a = pointX.g / u_R;",
                        "b = (pointZ.g - u_Z) / u_R;",

                        "normal = u_I * vec4(sign(b) * normalX.g, 0.0, normalZ.g, 0.0);",
                        "float green;",

                        "if(a > 2.0 || b > 2.0){",
                            "green = dot(normal, texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0)));",
                        "}else{",
                            "green = dot(normal,texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0)));",
                        "}",

                        "a = pointX.b / u_R;",
                        "b = (pointZ.b - u_Z) / u_R;",

                        "normal = u_I * vec4(sign(b) * normalX.b, 0.0, normalZ.b, 0.0);",
                        "float blue;",

                        "if(a > 2.0 || b > 2.0){",
                            "blue = dot(normal, texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0)));",
                        "}else{",
                            "blue = dot(normal,texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0)));",
                        "}",

                        "a = pointX.a / u_R;",
                        "b = (pointZ.a - u_Z) / u_R;",

                        "normal = u_I * vec4(sign(b) * normalX.a, 0.0, normalZ.a, 0.0);",
                        "float alpha;",

                        "if(a > 2.0 || b > 2.0){",
                            "alpha = dot(normal, texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0)));",
                        "}else{",
                            "alpha = dot(normal,texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0)));",
                        "}",

                        "gl_FragColor = vec4(red, green, blue, alpha);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_shape_half" : B_loop_half,
            "u_shape_tenth" : B_loop_tenth,
            "u_pointsX" : pointsX_tex_b,
            "u_pointsZ" : pointsZ_tex_b,
            "u_normalsX" : normalsX_tex_b,
            "u_normalsZ" : normalsZ_tex_b,
            "u_R" : spec.radius,
            "u_Z" : spec.height,
            "u_I" : -spec.current
        }).draw({
            triangles : 6,
            target : boundary_b
        }).set({
            "u_Z" : 0.0,
            "u_I" : spec.current
        }).draw({
            triangles : 6,
            target : boundary_b,
            blend : ['ONE', 'ONE']
        });

        // the matrix we will need to solve Ax = b
        // compute matrix A
        var boundary_A = webgl.addFrameBuffer({
            width: n_loops,
            height: n_loops,
            useFloat: true
        });

        // ---------------------------------------------------------------------
        // program for computing field normal contributions at points from variable loops
        var programVariableLoop = webgl.linkProgram({
            vertexShaderSource : render_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_shape_half;",
                    "uniform sampler2D u_shape_tenth;",

                    "uniform sampler2D u_loops;",
                    "uniform sampler2D u_points;",
                    "uniform sampler2D u_normals;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec4 loop = texture2D(u_loops, v_texCoord);",
                        "vec4 point = texture2D(u_points, v_texCoord);",
                        "vec4 normal = texture2D(u_normals, v_texCoord);",

                        // loop 0
                        "float a = point.x / loop.x;",
                        "float b = (point.z - loop.y) / loop.x;",

                        "vec4 tmp;",
                        "vec4 field = vec4(0.0);",

                        "if(a > 2.0 || abs(b) > 2.0){",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",

                        "field += tmp;",

                        // loop 0 reflected vertically
                        "a = point.x / loop.x;",
                        "b = (point.z - (1.0 - loop.y)) / loop.x;",

                        "if(a > 2.0 || abs(b) > 2.0){",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",

                        "field += tmp;",

                        // loop 1
                        "a = point.x / loop.z;",
                        "b = (point.z - loop.w) / loop.z;",


                        "if(a > 2.0 || abs(b) > 2.0){",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = - vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",

                        "field += tmp;",

                        // loop 1 reflected vertically
                        "a = point.x / loop.z;",
                        "b = (point.z - (1.0 - loop.w)) / loop.z;",


                        "if(a > 2.0 || abs(b) > 2.0){",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_tenth, vec2(a / 10.0, abs(b) / 10.0));",
                        "}else{",
                            "tmp = vec4(sign(b), 1.0, 1.0, 1.0) * texture2D(u_shape_half, vec2(a / 2.0, abs(b) / 2.0));",
                        "}",

                        "field += tmp;",

                        //  find normal to surface.
                        "gl_FragColor = vec4(dot(field, normal), 0.0, 0.0, 0.0);",

                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_shape_half" : B_loop_half,
            "u_shape_tenth" : B_loop_tenth,
            "u_loops" : loops_tex_A,
            "u_points" : points_tex_A,
            "u_normals" : normals_tex_b
        }).draw({
            triangles : 6,
            target : boundary_A
        });



        var result = eq.set_matrix(boundary_A).set_b(boundary_b).solve({
            tolerance : 1E-3,
            substep : 1,
            max_iterations : 10
        });

        var currents = result.result;

        // compute field everywhere using boundary solution
        for(loop = 0; loop < n_loops; loop++) {

            programCurrentLoop.set({
                "u_R0" : loops_arr_A[4 * loop], // loop0.x
                "u_Z0" : loops_arr_A[4 * loop + 1], // loop0.z
                "u_R1" : loops_arr_A[4 * loop + 2], // loop1.x
                "u_Z1" : loops_arr_A[4 * loop + 3], // loop1.z
                "u_I" : currents[loop]
            }).draw({
                triangles : 6,
                target : B,
                blend : ['ONE', 'ONE']
            });
        }

    };

    return exports;
});
