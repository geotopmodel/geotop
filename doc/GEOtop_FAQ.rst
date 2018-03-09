GEOtop FAQ
=============

:date: first version January 2017/ uopdate March 2018


Question: How can I restart my simulation ? 
--------------------------------------------

Answer: you need to set two keywords::

 SubfolderRecoveryFiles="rec"
 ContinuousRecovery=1

The first keyword sets the path to a folder (in this case called "rec") where GEOtop stores the recovery files. Make sure to create the directory if is not already present.
The second keywords sets the frequency in time when GEOtop saves the state variables in the recovery files. It is in "days", therefore if you want to save the state variables every 12 hours, you need to put 0.5, if you want to save them every day just put 1.
If the simulation suddenly breaks, you just need to re-launch it and it will start from the latest recovery time.
Please note that if you use the recovery files and you are re-starting a simulation, don't change the InitDateDDMMYYYYhhmm keyword, as the time variables are always set from that initial time.

[contributed by Matteo Dell’Amico]


Question: I got this error: ``[ERROR]: Point #   1 is out of the domain`` even if all my points are within the domain ! 
-----------------------------------------------------------------------------------------------------------------------
Answer: This generally indicated that METEOIO lib has not been compiled with ARC plugin enabled. IF this is the reason you should see something like that on your screen:: 

 [ERROR] MeteoIO: [IOHandler.cc:271] Cannot find plugin ARC as requested in file /XXXXX/geotop/tests/3D/panola/io_it.ini. Has it been    activated through ccmake? Is it declared in IOHandler::getPlugin?



 
 
