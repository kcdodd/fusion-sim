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

    exports.makeSORIterative = function(spec) {
        util.validate_object(spec, {
            n_power : 'number', // length(vector) = 4 * (2^n_power)^2,
            relaxation : [,'number'],
            webgl : [,'object']
        });

        var out = {};

        var vec_height = Math.pow(2, spec.n_power);

        var n_vec = vec_height * vec_height;

        var vec_length = 4 * n_vec;

        var mat_height = 2 * vec_height * vec_height;

        var n_mat = mat_height * mat_height;

        out.vec_length = vec_length;
        out.vec_height = vec_height;

        var omega = spec.relaxation || 1.0;

        // create canvas element for webgl to work on
        var webgl;

        if (spec.webgl) {
            webgl = spec.webgl;
        }else{
            var canvas = document.createElement("CANVAS");
            canvas.id = "webgl_canvas_JacobiIterative";
            canvas.width = vec_height;
            canvas.height = vec_height;
            canvas.style.display = "none";
            out.canvas = canvas;

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

        var x_set_arr = new Float32Array(vec_length);

        var x_set_tex = webgl.addTextureArray({
            width: vec_height,
            height: vec_height,
            array: x_set_arr,
            useFloat: true
        });

        var x_guess = webgl.addFrameBuffer({
            width : vec_height,
            height : vec_height,
            useFloat : true
        });

        var x_result = webgl.addFrameBuffer({
            width : vec_height,
            height : vec_height,
            useFloat : true
        });

        var x_stats = webgl.addFrameBuffer({
            width : vec_height,
            height : vec_height,
            useFloat : true
        });

        var b_set_arr = new Float32Array(vec_length);

        var b_set_tex = webgl.addTextureArray({
            width: vec_height,
            height: vec_height,
            array: b_set_arr,
            useFloat: true
        });

        var x_tex = null, A_tex = null, b_tex = null;

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
                    "precision highp float;",

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
        var programR = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",
                    // the original matrix where only red channel is used
                    // and rows extend the width of the texture, instead of being square
                    "uniform sampler2D u_A;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",

                        // pixel location from texture coordinates
                        "vec2 n = " + N(mat_height) + " * v_texCoord - vec2(0.5);",

                        // compute the original row and column from v_texCoord
                        "float row = floor(n.x / " + N(vec_height) + ") + " + N(2 * vec_height) + " * floor(n.y / " + N(vec_height) + ");",
                        "float col = 4.0 * (mod(n.x, " + N(vec_height) + ") + " + N(vec_height) + " * mod(n.y, " + N(vec_height) + "));",

                        // convert to texture coordinate \in [0,1] into the source matrix texture
                        "float c = " + N(1.0 / vec_length) + ";",

                        // the red component is the first element at the column
                        "float R = (row == col) ? 0.0 : -texture2D(u_A, c * vec2(col, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // green 2nd column
                        "float G = (row == col + 1.0) ? 0.0 : -texture2D(u_A, c * vec2(col + 1.0, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // blue 3rd column
                        "float B = (row == col + 2.0) ? 0.0 : -texture2D(u_A, c * vec2(col + 2.0, row)).r / texture2D(u_A, c * vec2(row, row)).r;",
                        // alpha 4th column
                        "float A = (row == col + 3.0) ? 0.0 : -texture2D(u_A, c * vec2(col + 3.0, row)).r / texture2D(u_A, c * vec2(row, row)).r;",

                        // pack columns into colors
                        "gl_FragColor = " + ((omega !== 1.0) ? N(omega) + " * " : "") + "vec4(R, G, B, A);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        // reconfigures the input vector b into the iteration vector C
        var programC = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",
                    // the original matrix where only red channel is used
                    // and rows extend the width of the texture, instead of being square
                    "uniform sampler2D u_A;",
                    "uniform sampler2D u_b;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        // get pixel number
                        "vec2 n = " + N(vec_height) + " * v_texCoord - vec2(0.5);",
                        // convert to row number
                        "float row = 4.0 * (n.x + " + N(vec_height) + " * n.y);",

                        "vec4 b = texture2D(u_b, v_texCoord);",
                        "float R = b.r / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row, row)).r;",
                        "float G = b.g / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 1.0, row + 1.0)).r;",
                        "float B = b.b / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 2.0, row + 2.0)).r;",
                        "float A = b.a / texture2D(u_A, " + N(1 / vec_length) + " * vec2(row + 3.0, row + 3.0)).r;",

                        "gl_FragColor = " + ((omega !== 1.0) ? N(omega) + " * " : "") + "vec4(R, G, B, A);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
        });

        // ---------------------------------------------------------------------
        // multiplies every element of a matrix by a component of a column vector
        var programMVproduct = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_M;",
                    "uniform sampler2D u_v;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        // compute the position in the vector from position in matrix
                        "vec2 row = " + N(2 * vec_height) + " * v_texCoord;",
                        "row = row - floor(row);",

                        // perform the matrix-vector product and stor individual results
                        "gl_FragColor = texture2D(u_M, v_texCoord) * texture2D(u_v, row);",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates,
        });

        // will store the raw result of the full product
        var mv_product = webgl.addFrameBuffer({
            width : mat_height,
            height : mat_height,
            useFloat : true
        });

        // reduces the number of x elements by 1/2 by summing adjacent elements
        // will be used successively to sum down to the resulting vector
        var sum_frag = function(num_x) {
            var src_arr = [
                "precision highp float;",
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

        for(var i = 0; i < spec.n_power; i++) {

            // ---------------------------------------------------------------------
            // program for summing adjacent elements
            //
            sum_programs[i] = webgl.linkProgram({
                vertexShaderSource : compute_vert(),
                fragmentShaderSource : sum_frag(mat_height / Math.pow(2, i))
            }).set({
                "a_position" : vertex_positions,
                "a_texCoord" : texture_coordinates,
                // the initial sum uses the raw result, while others use previous sums
                "u_M" : ((i > 0) ? sum_buffers[i-1] : mv_product)
            });

            // will hold the summation result, 1/2 the size of the previous texture in both directions
            sum_buffers[i] = webgl.addFrameBuffer({
                width : mat_height / Math.pow(2, i + 1),
                height : mat_height / Math.pow(2, i + 1),
                useFloat : true
            });
        }

        // ---------------------------------------------------------------------
        // finishes the sum and re-packs the resulting vector in the original form
        var programResult = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_Vsum;",
                    "uniform sampler2D u_C;",
                    "uniform sampler2D u_X;",

                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec2 offsetx = vec2(" + N(0.5/ (2 * vec_height)) + ", 0.0);",
                        "vec2 offsety = vec2(0.0, " + N(0.5/ (2 * vec_height)) + ");",

                        "float sumR = dot(texture2D(u_Vsum, v_texCoord - offsetx - offsety), vec4(1.0));",
                        "float sumG = dot(texture2D(u_Vsum, v_texCoord + offsetx - offsety), vec4(1.0));",
                        "float sumB = dot(texture2D(u_Vsum, v_texCoord - offsetx + offsety), vec4(1.0));",
                        "float sumA = dot(texture2D(u_Vsum, v_texCoord + offsetx + offsety), vec4(1.0));",

                        "gl_FragColor = vec4(sumR, sumG, sumB, sumA) + texture2D(u_C, v_texCoord) " + ((omega !== 1.0) ? (" + " + N(1.0 - omega) + " * texture2D(u_X, v_texCoord)") : "") + ";",
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
        // computes statistic between two textures
        var programStats = webgl.linkProgram({
            vertexShaderSource : compute_vert(),
            fragmentShaderSource : (function() {
                var src_arr = [
                    "precision highp float;",

                    "uniform sampler2D u_X1;",
                    "uniform sampler2D u_X2;",
                    // the texCoords passed in from the vertex shader.
                    "varying vec2 v_texCoord;",

                    "void main() {",
                        "vec4 x1 = texture2D(u_X1, v_texCoord);",
                        "vec4 x2 = texture2D(u_X2, v_texCoord);",
                        "vec4 diff = abs(x2 - x1);",
                        // the 0.25 is so the value does not exceed 1.0 for readback
                        // assuming x1 and x2 are [0,1]
                        "gl_FragColor = vec4(dot(x1, x2) * 0.25, dot(x1, x1) * 0.25, dot(x2, x2) * 0.25, max(max(max(diff.r, diff.g), diff.b), diff.a));",
                    "}"
                ];

                return src_arr.join('\n');
            })()
        }).set({
            "a_position" : vertex_positions,
            "a_texCoord" : texture_coordinates
        });

        // ---------------------------------------------------------------------
        // set matrix using a row-major 2D array or texture in red channel
        out.set_matrix = function(matrix) {

            if (Array.isArray(matrix)) {

                for(var row = 0; row < vec_length; row++) {
                    for(var col = 0; col < vec_length; col++) {
                        m_set_arr[4 * (col + row * vec_length)] = matrix[row][col];
                    }
                }

                m_set_tex.update();

                A_tex = m_set_tex;

            }else{
                // if it's not an array, try to use it as a texture.
                A_tex = matrix;
            }

            return out;
        };

        // ---------------------------------------------------------------------
        // set the b vector using an array or texture
        out.set_b = function(b) {

            if (Array.isArray(b)) {

                for(var j = 0; j < vec_length; j++) {
                    b_set_arr[j] = b[j];
                }

                b_set_tex.update();

                b_tex = b_set_tex;

            }else{
                b_tex = b;
            }

            return out;
        };

        // ---------------------------------------------------------------------
        // initalize the solution vector with an array or texture with components packed
        out.init_vector = function(vector) {

            if (Array.isArray(vector)) {

                for(var j = 0; j < vec_length; j++) {
                    x_set_arr[j] = vector[j];
                }

                x_set_tex.update();

                // render as the current solution result
                programSet.set({
                    'u_value' : x_set_tex
                }).draw({
                    triangles : 6,
                    target : x_result
                });



            }else{
                // try using vector as a texture and render directly to initial solution
                programSet.set({
                    'u_value' : vector
                }).draw({
                    triangles : 6,
                    target : x_result
                });
            }

            return out;
        };

        // ---------------------------------------------------------------------
        // compute a matrix-vector product, summing each row, and adding an offset.
        // x = R * x + C
        out.mv_product = function(target) {


            programMVproduct.draw({
                triangles : 6,
                target : mv_product
            });

            for(var i = 0; i < spec.n_power; i++) {
                sum_programs[i].draw({
                    triangles : 6,
                    target : sum_buffers[i]
                });
            }

            programResult.draw({
                triangles : 6,
                target : target
            });



            return out;
        };

        // ---------------------------------------------------------------------
        // iteratively solve

        var stats_arr = new Float32Array(vec_length);
        var x1_arr = new Float32Array(vec_length);
        var x2_arr = new Float32Array(vec_length);

        out.solve = function (params) {
            util.validate_object(params, {
                tolerance : 'number', // maximum relative difference between iterations
                substep : [,'number'], // take a certain number of iterations before computing statistics
                max_iterations : [, 'number'] // stop iterations if no convergence after this many steps
            });

            programSet.set({
                "u_value" : A_tex
            }).draw({
                triangles : 6
            });

            programR.set({
                "u_A" : A_tex
            }).draw({
                triangles : 6,
                target : R
            });

            programC.set({
                "u_A" : A_tex,
                "u_b" : b_tex
            }).draw({
                triangles : 6,
                target : C
            });

            // debug
            C.readPixels(x2_arr);
            console.log("C " + x2_arr);

            var m = new Float32Array(4 * n_mat);
            R.readPixels(m);
            console.log("R " + m);

            programMVproduct.set({
                "u_M" : R,
                "u_v" : x_guess
            });

            programResult.set({
                "u_C" : C,
            });

            if (omega !== 1.0) {
                programResult.set({
                    "u_X" : x_guess,
                });
            }

            programStats.set({
                "u_X1" : x_guess,
                "u_X2" : x_result
            });

            programSet.set({
                'u_value' : x_result
            });



            var correlation = 0.0;
            var max_diff = 0.0;


            var x1 = 0;
            var x2 = 0;
            var x1x2 = 0;
            var x1x1 = 0;
            var x2x2 = 0;
            var diff = params.tolerance + 1;
            var iteration = 0;

            // keep iterating until the guess and the result are within tolerance
            while(iteration < params.max_iterations && diff > params.tolerance) {


                for(var sub = 0; sub < (params.substep || 1); sub++) {
                    // update the guess
                    programSet.draw({
                        triangles : 6,
                        target : x_guess
                    });

                    // copmute new result
                    out.mv_product(x_result);

                }


                // compute statistics
                programStats.draw({
                    triangles : 6,
                    target : x_stats
                });

                x_stats.readPixels(stats_arr);
                x_guess.readPixels(x1_arr);
                x_result.readPixels(x2_arr);

                max_diff = 0.0;

                for(i = 0; i < n_vec; i++) {
                    x1 += x1_arr[4 * i] + x1_arr[4 * i + 1] + x1_arr[4 * i + 2] + x1_arr[4 * i + 3];
                    x2 += x2_arr[4 * i] + x2_arr[4 * i + 1] + x2_arr[4 * i + 2] + x2_arr[4 * i + 3];
                    x1x2 += stats_arr[4 * i];
                    x1x1 += stats_arr[4 * i + 1];
                    x2x2 += stats_arr[4 * i + 2];

                    max_diff = Math.max(max_diff, stats_arr[4 * i + 3]);
                }

                // Pearson product-moment correlation coefficient
                correlation = (vec_length * x1x2 - x1 * x2) / Math.sqrt((vec_length * x1x1 - x1 * x1) * (vec_length * x2x2 - x2 * x2));
                diff = 2 * vec_length * max_diff / (Math.abs(x1) + Math.abs(x2));

                console.log(iteration + ":" + x2_arr);

                iteration++;
            }

            return {
                correlation : correlation,
                diff : diff,
                iterations : iteration,
                result : x2_arr
            };

        };

        // ---------------------------------------------------------------------
        out.x_result_tex = function() {

            return x_result;
        };



        return out;
    };



    return exports;
});
