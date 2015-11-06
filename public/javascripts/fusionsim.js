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



define([
    'angular',
    'utilities',
    'empic',
    'matrix_webgl'
], function (angular, utilities, empic, matrix_webgl){
    "use strict";

    var app = angular.module('fusionsim', []);

    app.controller('simulation', [
        '$scope',
        '$interval',
        function($scope, $interval){

/*
            var eq = matrix_webgl.makeSORIterative({
                n_power : 1,
                relaxation : 1.0
            });

            var A = [];
            var b = [];
            var x = [];

            for(var row = 0; row < 16; row++) {
                A[row] = [];

                for(var col = 0; col < 16; col++) {
                    A[row][col] = 0.0;
                }

                A[row][row] = Math.random();
                b[row] = Math.random();
            }

            eq.set_matrix(A).set_b(b);

            var result = eq.solve({
                tolerance : 1E-3,
                substep : 1,
                max_iterations : 100
            });

            console.log("done: " + result.result);

            return;
*/

            $scope.run = false;
            $scope.fps = 0;

            var nparticles = 160000;

            var spec = {
                radius : 1,
                height : 2,
                nr : 400,
                nz : 800,
                dt : 2E-9,
                nparticles : 400,
                particle_mass : 1.67e-27, //kg : proton
                particle_charge : 1.602e-19 // C
            };

            var simulation = empic.makeCylindricalParticlePusher(spec);

            var init_position = [];
            var init_velocity = [];
            var sink = [];
            var source = [];

            var i, j, k;

            for(i = 0; i < spec.nr; i++) {
                source[i] = [];
                sink[i] = [];

                for(j = 0; j < spec.nz; j++) {
                    sink[i][j] = 1.0;
                    source[i][j] = 0.0;
                }
            }

            // sink mask
            for(j = 0; j < spec.nz; j++) {
                //sink[0][j] = 0;
                sink[spec.nr-1][j] = 0;
            }
            for(i = 1; i < spec.nr-1; i++) {
                sink[i][0] = 0;
                sink[i][spec.nz-1] = 0;
            }

            // source pdf

            for(i = 0; i < 50; i++) {

                for(j = 350; j < 450; j++) {

                    source[i][j] = 1.0;
                }
            }

            // particles
            for(i = 0; i < nparticles; i++) {
                init_position[i] = [0.2 * (Math.random() - 0.5), 0.2 * (Math.random() - 0.5), 0.2 * (Math.random() - 0.5) + 1];
                init_velocity[i] = [0.002 * (Math.random()-0.5), 0.002 * (Math.random()-0.5), 0.002 * (Math.random()-0.5)];
            }

            simulation.set({
                position : init_position,
                velocity : init_velocity,
                sink_mask : sink,
                source_pdf : source
            });

            simulation.addCurrentLoop(0.8, 2.0, -10000000);
            simulation.addCurrentLoop(0.8, 0.0, 10000000);

            //simulation.addCurrentLoop(0.5, 1.0, 10000000);
            //simulation.addCurrentZ(5000000);
            //simulation.addBZ(0.01);

            //simulation.addBTheta(0.01);

            //simulation.addSpindleCuspPlasmaField(1.0, 0.5);

            simulation.precalc();

            var plot = document.getElementById("plot");
            var plot_ctx = plot.getContext("2d");

            simulation.density();
            plot_ctx.drawImage(simulation.canvas, 0, 0);

            console.log("initialized");

            window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
                    window.webkitRequestAnimationFrame || window.oRequestAnimationFrame;


            $scope.start = function() {
                console.log("started");

                var t_last = Date.now();
                var N_frames = 0;

                $scope.run = true;

                var step = function(){

                    simulation.step();

                    simulation.density();

                    plot_ctx.fillStyle = "rgba(0,0,0,0)";
                    plot_ctx.clearRect(0,0,spec.nr, spec.nz);
                    plot_ctx.drawImage(simulation.canvas, 0, 0);

                    if ($scope.run) {

                        var t_cur = Date.now();
                        N_frames++;

                        if (t_cur - t_last > 1000){
                            $scope.$apply(function(){

                                $scope.fps = (N_frames * 1000/(t_cur - t_last)).toFixed(0);
                                t_last = t_cur;
                                N_frames = 0;
                            });
                        }


                        requestAnimationFrame(step);
                    }else{
                        $scope.$apply(function(){
                            $scope.fps = 0;
                        });

                    }
                };

                requestAnimationFrame(step);
            };

            $scope.stop = function() {
                console.log("stopped");
                $scope.run = false;
            };
	}]);

    app.controller('MainCtrl', [
        '$scope',
        '$interval',
        function($scope, $interval){

	}]);

});
