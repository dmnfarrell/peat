## DataPipeline for programmers ##

DataPipeline is written in the python language and can be extended by writing a custom Importer class or adding filters/functions. Users who know [Python](http://www.python.org/) and are interested in adding their own importers may find the information here useful. The code is compatible with Python versions>=2.5.

## Classes ##

The (somewhat non-standard) diagram below shows how the core classes relate to the main components. We currently import several plotting and fitting methods from PEATDB, so this package is a dependency when installing the source package. Note that PEATDB has itself other dependencies such as ZODB, but these are not required for DataPipeline to work.

<img src='http://peat.googlecode.com/svn/wiki/images/datapipeline_classes.png' width='800'>

The core class is <code>Base.Pipeline</code> which handles all file reading/writing, queuing, fitting and joining of data, managing configuration files.<br>
<br>
The <code>BaseImporter</code> class does the actual import from the file, given a configuration. It has basic methods for getting rows, headers etc. This is always subclassed and the doImport method overridden to achieve the desired functionality. There are currently 8 Importer classes built in to DataPipeline that handle the most common text formats.<br>
These are:<br>
<br>
<ul><li><code>DatabyColImporter</code>
</li><li><code>DatabyRowImporter</code>
</li><li><code>PairedDatabyColImporter</code>
</li><li><code>PairedDatabyRowImporter</code>
</li><li><code>GroupedDatabyColImporter</code>
</li><li><code>GroupedDatabyRowImporter</code></li></ul>

See DataPipelineFormats<br>
<br>
<h2>Adding a custom importer</h2>

To add your own custom importers, simple subclass the <code>BaseImporter</code> class and put the code in Custom.py so that the importer will be loaded automatically when required.<br>
<br>
<ul><li>Set the 'name' property of the importer, this label will be used in the <i>format</i> keyword in the configuration file to identify your importer.<br>
</li><li>Override the doImport method at minimum.<br>
</li><li>Look at the other basic importers in Importer.py to see how the various methods are used such as getRowHeader() and so on. You can override these other methods if required, but generally they should be left alone.</li></ul>

A template for a new importer would look like this:<br>
<br>
<pre><code>class MyCustomImporter(BaseImporter):<br>
    """This importer handles data formatted in this format.. """<br>
<br>
    name = 'mycustomimporter'<br>
<br>
    def __init__(self, cp):<br>
        BaseImporter.__init__(self, cp)<br>
        return<br>
<br>
    def doImport(self, lines):<br>
        #your code here<br>
</code></pre>

<h2>Adding filters/functions</h2>

Pre-defined filters are implemented by the <code>Processor</code> class. The <code>doProcessingStep()</code> method of the Pipeline class (called when the <code>run()</code> method is executed) first checks that the names are found in the <code>predefined</code> attribute of the Processor class, if not the method returns.<br>
<br>
<pre><code>X = Processor(self.functionsconf)<br>
names = [i[1] for i in self.functions]<br>
for n in names:<br>
    if n not in X.predefined:<br>
        print 'function %s not found' %n<br>
        return data<br>
data = X.doFunctions(names, data)<br>
</code></pre>

<code>Processor.doFunctions</code> then executes the filters simply by looping over the function names in the order they are defined in the configuration file.<br>
<br>
For now, functions can simply be added by inserting new methods into the Processor class found in Processing.py. Also add the name of the function to the classes predefined attribute. New functions should be specified in the following form:<br>
<br>
<pre><code>def somefunction(self, x, y):<br>
    #process values here, if x or y is unchanged, just return it<br>
    return newx,newy<br>
</code></pre>