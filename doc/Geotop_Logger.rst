GEOtop Logger and GlobalLogger usage
====================================


:author: Gianfranco Gallizia (eXact-lab srl)
:date: september 2014

GEOtop logging facility is based on two classes:

* ``Logger``: a multi-stream, severity-based logger.
* ``GlobalLogger``: a singleton that handles the main Logger.

To use the ``GlobalLogger`` you have to get the instance and then you can
start logging::

    #include "global_logger.h"

    void myFunc() {
        geotop::logger::GlobalLogger* lg =
            geotop::logger::GlobalLogger::getInstance();

        lg->log("Message");
    }

This will prompt on the console and on the main log file the following
message::

    [NOTICE]: Message

It is possible to select the message's severity among the following:

* ``TRACE`` for maximum verbosity.
* ``DEBUG`` for debug messages.
* ``NOTICE`` standard severity.
* ``WARNING`` e.g. when GEOtop uses default values instead of user-provided values.
* ``ERROR`` strong but (mostly) recoverable errors.
* ``CRITICAL`` unrecoverable errors.

All these levels are defined in ``logger.h`` as an ``enum`` called
``geotop::logger::severity_levels``. To select the severity level use the
``log`` method with one of the ``severity_levels`` as second parameter::

    lg->log("Message", geotop::logger::WARNING);

The above call will output this::

    [WARNING]: Message

Formatted output
^^^^^^^^^^^^^^^^

There are two methods to generate ``printf``-style formatted messages:

* ``logf`` which uses the default severity (NOTICE).
* ``logsf`` which uses the provided severity.

These are their signatures::

    void logf(const char* format, ...);
    void logsf(severity_levels severity, const char* format, ...);

E.G.::

    lg->logf("Testing: %d", 1);
    lg->logsf(geotop::logger::WARNING, "Testing: %d", 2);

Will output::

    [NOTICE] Testing: 1
    [WARNING] Testing: 2

Verbatim output
^^^^^^^^^^^^^^^

You can print a message without pre-pending a severity with the following
methods:

* ``writeAll`` writes a ``std::string`` on all the output streams managed by the logger.
* ``writefAll`` acts like a ``printf`` on all the output streams managed by the logger.


Status of the implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the moment  meteo.cc meteodistr.cc and input.cc are using the logger. Please refer to them how to use it. 
