Minor-Gravity Package
===============================
A package to compute irregular gravitational fields around minor bodies, such as Asteroids and Comets.
---------------------------------

I've implemented a modified version of the original POLYHEDRON code from D. Tsoulis into the Minor-Gravity package, that can compute the gravitational potential, and its first and second derivatives of a homogenous polyhedron, according to Petrovic (J of G, 1996). In addition, the Minor-Gravity package also can compute gravitational features of a polyhedron using the mascons approach.

I've made some adaptations in the original Tsoulis (2012) code. For example, Minor-Gravity package writes the output file in a compressed format file using ASCII characters through subroutine ``mio_out.for``, adapted from Mercury package (Chambers, 1999). The compressed output file can have its data decompressed in a human-readable file and the user can also choose slices of boxes with arbitrary dimensions (for example, projection planes xOy, xOz and yOz) from original huge input grid. This can be made running the ``elementgr.for`` source from Minor-Gravity package. This reduces the original output file approximately ``3 times less`` considering double precision floating numbers (high output precision). For example, a total simulation size for KBO Ultima Thule is about ``6.7 GB``, while without compressed format file it will have an approximated ``19.5 GB`` size, specially using original ``polyhedron.f`` code. Another advantage is that Minor-Gravity package can compute numerically gravitational features from an irregular polyhedral shape at a bunch of choosen points, differently of the original polyhedron.f code that computes the gravity field always at center of mass of the polyhedron. This makes the algorithm more useful to deal numerically with this kind of management.

Notable contents of this repository
---------------------------

*    ``minor-gravity.for:`` The main program file.  It requires ``polyhedron.inc`` to compile.
*    ``command.for:`` Additional source file to get overall CPU % usage.
*    ``dhms.for:`` Additional source file to find the approximately remaining time of the integration.
*    ``get_pid.for:`` Additional source file to get PID process of the main program.
*    ``timestamp2.for:`` Additional source file to print out the current YMDHMS date as a timestamp.
*    ``./output/elementgr.for:`` This programme converts an output file created by main program file into a set of file(s) containing the gravitational features of the polyhedron. These file(s) can be used as the basis for making graphs or movies using a graphics package as GNUPlot.

 Input Files
 ---------------

*    ``files.in:`` This should contain a list of names for the input and output files used by main program file. List each file name on a separate line.
*    ``polyhedron.in:`` This file contains parameters that control how an integration is carried out. Any lines beginning with the ) character are assumed to be comments, and will be ignored by main program file. All initial parameters of this file are set for the application of the KBO Arrokoth.
*    ``./in/vertex.in:`` This file contains the coordinates x, y and z of the vertices of the polyhedron.
*    ``./in/face.in:`` This file contains the index of the vertices belonging at each given facet, in such an order that the plane normal is always pointing outside the body.
*    ``./in/cube.in:`` (optional) This file is used from gravitational potential of a cubic grid computed previously through the mascons technique.
*    ``./in/grid.in:`` (optional) This file is used from a cubic grid generated previously.
*    ``./in/message.in:`` This file contains the text of various messages output by main program file, together with an index number and the number of characters in the string (including spaces used for alignment).
*    ``./output/element.in:`` This file contains parameters and options used by elementgr.for. Any lines beginning with the ) character are assumed to be comments, and will be ignored by the source.

 Dump files
 ---------------

*    ``polyhedron.dmp:`` A dump file containing the integration parameters. If your computer crashes while main program file is doing an integration, all is not lost. You can continue the integration from the point at which main program file last saved data in dump files, rather than having to redo the whole calculation.
*    ``polyhedron.tmp:`` The same as above (backup file). If for some reason one of the dump files has become corrupted, look to see if you still have a set of files with the extension .tmp produced during the original integration (if you have subsequently used main program file to do another integration in the same directory, you will have lost these unfortunately). These .tmp files are duplicate copies of the dump files. Copy each one so that they form a set of uncorrupted dump files, and then run the compiled version of main program file.
*    ``restart.dmp:`` An additional dump file containing other variables used by main program.
*    ``restart.tmp:`` The same as above (backup file).

 N.B. It is important that you replace all the dump files with the .tmp files in this way, rather than just the file that is corrupted.

 Output files
 ---------------

*    ``./out/v1.out:`` Gravitational features of the polyhedron in the integration, produced at incremented intervals.

How to compile and run
----------------------

In the current directory, you should have all files needed to compile and execute the program. Don't hesitate to contact me if any problem arises when trying to use it.
Current package is preconfigured for Linux and Mac OS X users with compile and clean tasks. Use your favorite FORTRAN compiler, such as ``gfortran``, ``f77`` or ``ifort``, to create an executable.  For instance, on Linux or Mac, try:

   ``gfortran minor-gravity.for command.for dhms.for get_pid.for timestamp2.for -o minor-gravity``

There will likely be warnings due to the code being written in FORTRAN77, but it should compile.  Copy or link the executable wherever you want (wherever your input files are) to run your code using ``./minor-gravity``.

Or use ifort compiler for no warnings:

   ``ifort minor-gravity.for command.for dhms.for get_pid.for timestamp2.for -o minor-gravity``

Or use the ``Makefile:``

   To compile on Linux or Mac OS X just run ``make`` on the terminal at current directory. You need to have gfortran compiler. If you don't have you need to install it. For other compilers change the corresponding lines in ``Makefile``. This file can be changed if you want to use another compiler as the default one (gfortran). Usually I use the compiler ifort, which leads to faster integration time.

   To clean the directory *.o files run ``make clean``.

Or use the executable files ``minor-gravity_gfort`` or ``minor-gravity_ifort`` for Linux and Mac OS X compilers ``gfortran`` or ``ifort``, respectively:

   First type ``chmod +x minor-gravity_gfort`` and after that to run your code using ``./minor-gravity_gfort``.

Tricks and Caveats
------------------

Unfortunately, the code needs to be recompiled any time parameters in the ``polyhedron.inc`` file get changed.

Disclaimers
------------

* The changes have been successful tested with ``polyhedra`` and ``mascons`` approaches.
* I've fixed all the errors I've found.  If you find a bug, let me know so we can try to fix it.
* Any feedback is appreciated, especially bugs, suggestions, or possible contributions.
* Are you going to publish? Please acknowledge the use of my code in any publication referencing:

   ``Amarante, A., Winter, O. C. (2020), Surface Dynamics, Equilibrium Points and Individual Lobes of the Kuiper Belt Object (486958) Arrokoth. MNRAS.``

* The source codes of this repository are available on reasonable request.

