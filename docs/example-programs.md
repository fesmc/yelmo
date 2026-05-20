# Example programs

The Yelmo base code provides a static library interface that can be used in
other programs, as well as a couple of stand-alone programs for running
certain benchmarks. Here we provide more examples of how to use
Yelmo:

1. Program template to connect with other models/components.
2. Stand-alone ice sheet with full boundary forcing.

In both cases, it is necessary to download the Yelmo repository
separately, as well as compile the Yelmo static library
(see [Getting started](getting-started.md)).

## Program template

This is a minimalistic setup that allows you to run Yelmo with
no dependencies and a straightforward Makefile. This template
can be used to design a new stand-alone Yelmo experiment, or
to provide guidance when adding Yelmo to another program.

Clone the repository from [https://github.com/palma-ice/yelmot](https://github.com/palma-ice/yelmot)

## Stand-alone ice sheet with full boundary forcing (yelmox)

This setup is suitable for glacial-cycle simulations, future simulations
or any other typical (realistic) ice-sheet model simulation.

Clone the repository from [https://github.com/fesmc/yelmox](https://github.com/fesmc/yelmox)
