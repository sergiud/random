.. -*- restructuredtext -*-

What is random?
===============

random is a collection of various C++ pseudo-random number generators and
distributions that are intended to accompany the Boost Random Number Library.

Getting started
===============

random requires the Boost C++ Libraries 1.43. Previous versions might work as
well. To compile the tests make sure first CMake 2.8 is installed, then enter ::

  $ cmake . -DCMAKE_BUILD_TYPE=Release

in your terminal or command prompt on Windows inside project's directory to
generate the appropriate configuration that can be used to compile the tests
using make/nmake or inside an IDE.
