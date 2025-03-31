# DynamO:- A general event driven simulator.

## C++ Source Code

The source code is divided into several seperate projects, each with its own
directory.

`dynamo/`

This is where all the Event-Driven simulation specific code is and is probably
what you are interested in if you're reading this.

`magnet/`

Magnet is a header-only library, filled with useful and generic C++ functions
and classes. DynamO uses the magnet library for stuff like sorting, collision
detection, and input/output in XML.

`coil/`

Coil is a visualization library. It uses OpenCL and OpenGL to perform fast
real-time renderings of simulations. In DynamO, the coil library provides the
interactive visualiser.
