@echo on
::
::
:: Set Environmental Variables for Geant4.10.7.2
:: (this step is not needed if you have installed Geant from scratch)
::
call SetEnvironment.bat
::
::
:: Run TrueBq01 user program
::
:: Can also be run it from command prompt by typing "RelwithDebInfo\TrueBq01.exe batch.mac" (assuming Environment has been set)
:: Interactive mode (visualization) "RelwithDebInfo\TrueBq01.exe" (runs vis.mac, but you don't type it)
::
@cmd.exe /K ".\RelWithDebInfo\TrueBq01.exe batch.mac && echo[ && echo Completed RelWithDebInfo\TrueBq01.exe batch.mac && echo[ && echo Output in .out files && echo["
