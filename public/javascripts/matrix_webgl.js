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
/*
    matrix_webgl

    This module provides methods to help solve a linear system of equations using
    iteration: Ax = b

    computations are performed on a gpu using webgl.
*/
define(['utilities'], function (util){
    "use strict";

    var exports = {};

    // formatting numbers for shaders
    var N = function(num) {
        return num.toFixed(20);
    };

    exports.makeJacobiIterative = function(spec) {
        util.validate_object(spec, {
            n_power : 'number', // length(vector) = 4 * (2^n_power)^2
            webgl : [,'object']
        });

        var out = {};

        var vec_height = Math.pow(2, spec.n_power);

        var n_vec = vec_height * vec_height;

        var vec_length = 4 * n_vec;

        var mat_height = Math.pow(2, 2 * spec.n_power + 1);

        var n_mat = mat_height * mat_height;

        out.vec_length = vec_length;

        // create canvas element for webgl to work on
        var webgl;

        if (spec.webgl) {
            webgl = spec.webgl;
        }else{
            canvas = document.createElement("CANVAS");
            canvas.id = "webgl_canvas_CylindricalParticlePusher";
            canvas.width = 100;
            canvas.height = 100;
            canvas.style.display = "none";

            document.body.appendChild(canvas);

            webgl = util.webGL(canvas);
        }

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


        var m_set_arr = new Float32Array(4 * 4 * n_mat);

        var m_set_tex = webgl.addTextureArray({
            width: 2 * mat_height,
            height: 2 * mat_height,
            array: m_set_arr,
            useFloat: true
        });

        var x_set_arr = new Float32Array(4 * n_vec);

        var x_set_tex = webgl.addTextureArray({
            width: vec_height,
            height: vec_height,
            array: x_set_arr,
            useFloat: true
        });

        var x_result = webgl.addFrameBuffer({
            width : vec_height,
            height : vec_height,
            useFloat : true
        });

        var b_set_arr = new Float32Array(4 * n_vec);

        var b_set_tex = webgl.addTextureArray({
            width: vec_height,
            height: vec_height,
            array: c_set_arr,
            useFloat: true
        });

        var x_tex, A_tex, b_tex;

        // the iteration matrix
        var R = webgl.addFrameBuffer({
            width : mat_height,
            height : mat_height,
            useFloat : true
        });

        // the constant vector
        var C = webgl.addFrameBuffer({
            width : vec_height,
            height : vec_height,
            useFloat : true
        });

        //
        // General shader to compute a texture frame buffer.
        //
        var compute_vert = function() {
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
        }; // compute_vert()

        // ---------------------------------------------------------------------
        // program for setting value of a texture
        //
        var programSet = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
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
        // reconfigures the input matrix into the iteration matrix R
        var programRJacobi = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",
                    // the original matrix where only red channel is used
                    // and rows extend the width of the texture, instead of being square
                    "uniform sampler2D u_A;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        // compute the original row and column from v_texCoord
                        "vec2 n = " + N(mat_height) + " * v_texCoord;",
                        "float row = floor(n.x / " + N(vec_height) + ") + " + N(2 * vec_height) + " * floor(n.y / " + N(vec_height) + ");",
                        "float col = 4.0 * (mod(n.x, " + N(vec_height) + ") +  " + N(vec_height) + " * mod(n.y, " + N(vec_height) + "));",

                        // convert to texture coordinate \in [0,1] into the source matrix texture
                        "float c = " + N(1 / vec_length) + ";",

                        // the red component is the first element at the column
                        "float R = (row == col) ? vec4(0.0) : -texture2D(u_A, c * vec2(col, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // green 2nd column
                        "float G = (row == col + 1) ? vec4(0.0) : -texture2D(u_A, c * vec2(col + 1, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // blue 3rd column
                        "float B = (row == col + 2) ? vec4(0.0) : -texture2D(u_A, c * vec2(col + 2, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // alpha 4th column
                        "float A = (row == col + 3) ? vec4(0.0) : -texture2D(u_A, c * vec2(col + 3, row)).r / texture2D(u_A, c * vec2(row, row)).r;",

                        // pack columns into colors
                        "gl_FragColor = vec4(R, G, B, A);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_A" : A_tex
        });

        // ---------------------------------------------------------------------
        // reconfigures the input vector b into the iteration vector C
        var programCJacobi = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",
                    // the original matrix where only red channel is used
                    // and rows extend the width of the texture, instead of being square
                    "uniform sampler2D u_A;",
                    "uniform sampler2D u_b;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "float row = " + N(4 * vec_height) + " * (v_texCoord.x + " + N(vec_height) + " * v_texCoord.y);",
                        "float R = texture2D(u_b, v_texCoord).r / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row, row)).r;",
                        "float G = texture2D(u_b, v_texCoord).g / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 1, row + 1)).r;",
                        "float B = texture2D(u_b, v_texCoord).b / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 2, row + 2)).r;",
                        "float A = texture2D(u_b, v_texCoord).a / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 3, row + 3)).r;",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_A" : A_tex,
            "u_b" : b_tex
        });

        // ---------------------------------------------------------------------
        // multiplies every element of a matrix by a component of a column vector
        var programMVproduct = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_M;",
                    "uniform sampler2D u_v;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec2 row = " + N(2 * vec_height) + " * v_texCoord;",
                        "row = row - floor(row);",
                        "gl_FragColor = texture2D(u_M, v_texCoord) * texture2D(u_v, row);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
        });

        var mv_product = webgl.addFrameBuffer({
            width : mat_height,
            height : mat_height,
            useFloat : true
        });

        // reduces the number of x elements by 1/2 by summing adjacent elements
        var sum_frag = function(num_x) {
            var src_arr = [
                "precision mediump float;",
                // matrix to be summed over
                "uniform sampler2D u_M;",
                // the texCoords passed in from the vertex shader.
                "varying vec2 v_texCoord;",

                "void main() {",
                    "vec2 offsetx = vec2(" + N(0.5/num_x) + ", 0.0);",
                    "vec2 offsety = vec2(0.0, " + N(0.5/num_x) + ");",
                    "gl_FragColor = texture2D(u_M, v_texCoord + offsetx + offsety) + texture2D(u_M, v_texCoord - offsetx + offsety) + texture2D(u_M, v_texCoord + offsetx - offsety) + texture2D(u_M, v_texCoord - offsetx - offsety);",
                "}"
            ];

            return src_arr.join('\n');
        }; // sum_frag

        var sum_buffers = [];
        var sum_programs = [];

        for(i = 0; i < spec.n_power; i++) {

            // ---------------------------------------------------------------------
            // program for summing adjacent elements
            //
            sum_programs[i] = webgl.linkProgram({
                vertexShaderSource : compute_vert(),
                fragmentShaderSource : sum_frag(Math.pow(2, 2 * spec.n_power + 1 - i))
            }).set({
                "a_position" : vertex_positions,
                "a_texCoord" : texture_coordinates,
                "u_M" : ((i > 0) ? sum_buffers[i-1] : mv_product)
            });

            sum_buffers[i] = webgl.addFrameBuffer({
                width : Math.pow(2, 2 * spec.n_power - i),
                height : Math.pow(2, 2 * spec.n_power - i),
                useFloat : true
            });
        }

        // ---------------------------------------------------------------------
        // finishes the sum and re-packs the resulting vector in the original form
        var programResult = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision mediump float;",

                    "uniform sampler2D u_Vsum;",
                    "uniform sampler2D u_C;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec2 offsetx = vec2(" + N(0.5/Math.pow(2, 2 * spec.n_power + 1)) + ", 0.0);",
                        "vec2 offsety = vec2(0.0, " + N(0.5/Math.pow(2, 2 * spec.n_power + 1)) + ");",

                        "float sumR = dot(texture2D(u_Vsum, v_texCoord - offsetx - offsety), vec4(1.0));",
                        "float sumG = dot(texture2D(u_Vsum, v_texCoord + offsetx - offsety)), vec4(1.0));",
                        "float sumB = dot(texture2D(u_Vsum, v_texCoord - offsetx + offsety), vec4(1.0));",
                        "float sumA = dot(texture2D(u_Vsum, v_texCoord + offsetx + offsety), vec4(1.0));",

                        "gl_FragColor = vec4(sumR, sumG, sumB, sumA) + texture2D(u_C, v_texCoord);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
            "u_Vsum" : sum_buffers[spec.n_power-1]
        });



        // ---------------------------------------------------------------------
        out.set_matrix = function(matrix) {
            for(k = 0; k < 4 * n_mat; k++) {

                var row = Math.floor(k / vec_length);
                var col = k % vec_length;

                m_set_arr[4 * k] = matrix[row][col];
            }

            m_set_tex.update();

            A_tex = m_set_tex;
        };

        // ---------------------------------------------------------------------
        out.set_matrix_tex = function(matrix_tex) {

            A_tex = matrix_tex;

        };

        // ---------------------------------------------------------------------
        out.set_b = function(b) {

            for(j = 0; j < n_vec; j++) {
                b_set_arr[4 * j] = b[4 * j];
                b_set_arr[4 * j + 1] = b[4 * j + 1];
                b_set_arr[4 * j + 2] = b[4 * j + 2];
                b_set_arr[4 * j + 3] = b[4 * j + 3];
            }

            b_set_tex.update();

            b_tex = b_set_tex;
        };

        // ---------------------------------------------------------------------
        out.set_b_tex = function(b) {
            b_tex = b;
        };

        // ---------------------------------------------------------------------
        out.init_vector = function(vector) {

            for(j = 0; j < n_vec; j++) {
                x_set_arr[4 * j] = vector[4 * j];
                x_set_arr[4 * j + 1] = vector[4 * j + 1];
                x_set_arr[4 * j + 2] = vector[4 * j + 2];
                x_set_arr[4 * j + 3] = vector[4 * j + 3];
            }
            x_set_tex.update();

            programSet.set({
                'u_value' : x_set_tex
            }).draw({
                triangles : 6,
                target : x_result
            });

        };

        // ---------------------------------------------------------------------
        out.init_vector_tex = function(vector_tex) {

            programSet.set({
                'u_value' : vector_tex
            }).draw({
                triangles : 6,
                target : x_result
            });

        };

        // ---------------------------------------------------------------------
        out.mv_product = function(target) {
            programMVproduct.draw({
                triangles : 6,
                target : mv_product
            });

            for(i = 0; i < num_sums; i++) {
                sum_programs[i].draw({
                    triangles : 6,
                    target : sum_buffers[i]
                });
            }

            programResult.draw({
                triangles : 6,
                target : target
            });
        };

        // ---------------------------------------------------------------------
        out.solve = function () {
            programMVproduct.set({

            });
        };

        // ---------------------------------------------------------------------
        out.x_result = function() {

            x_result.readPixels(x_set_arr);

            var result = [];

            for(j = 0; j < n_vec; j++) {
                result[4*j] = x_set_arr[4*j];
                result[4*j + 1] = x_set_arr[4*j + 1];
                result[4*j + 2] = x_set_arr[4*j + 2];
                result[4*j + 3] = x_set_arr[4*j + 3];
            }

            return result;
        };

        // ---------------------------------------------------------------------
        out.x_result_tex = function() {

            return x_result;
        };



        return out;
    };

    return exports;
});
