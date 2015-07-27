## About ##

PEATDB supports a simple plugin system for the client application that allows a python developer to add functionality. This may simply be some custom set of widgets for showing items in the database or a dialog for communicating with other software. Plugins should inherit from the Plugin class. A template source file is provided and the appropriate methods are overriden to add the required GUI elements. The PEAT API contains table widgets that allow rows from the database (or any arbitrary data) to be displayed in the plugin window. Plotting/fitting functionality can be added using the Ekin classes or by directly utilising the matplotlib library. Plugins source .py files are placed inside a plugin folder. These classes dynamically loaded at startup and shown in the main menu. Plugins can be reloaded whilst the application is running.

## Current Plugins ##

We have developed several plugins for our own analysis and make them available with PEAT:
  * PEATSA Plugin: Submit and retrieve structure-based energy calculations using PEATSA. See [PEATSAPlugin](PEATSAPlugin.md)
  * Correlation Analysis: A plugin we are developing that does analysis of experimental vs model prediction correlations. Requires scipy.
  * Van't Hoff Analysis: Do analysis on thermal denauration (CD) data using various methods from the literature.

See also [PEATDB](PEATDB.md)