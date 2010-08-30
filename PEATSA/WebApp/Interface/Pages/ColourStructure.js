/**
 Protein Engineering Analysis Tool Structure Analysis (PEATSA)
 Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
 
 Author: Michael Johnston
 
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
 
 Contact information:
 Email: Jens.Nielsen_at_gmail.com
 Normal mail:
 Jens Nielsen
 SBBS, Conway Institute
 University College Dublin
 Dublin 4, Ireland
 **/

//Functions for coloring the pdb structure displayed in Jmol on the results page
//It requires that JSGetCalculationData.php is included.

//Returns an array containing the data for calculation.
//If calculation is not finished this function returns nil
function dataForCalculation(calculation)
{
	if(calculation == "Stability")
	{
		return create_Stability_data_array();
	}
	else if(calculation == "Scan")
	{
		return create_Scan_data_array();
	}
	else if(calculation == "Binding")
	{
		return create_Binding_data_array();
	}
}

Array.prototype.max = function() 
{
	var max = this[0];
	var len = this.length;
	for (var i = 1; i < len; i++)
	{
		if (this[i] > max) 
		{
			max = this[i];
		}
	} 
	return max;
}

Array.prototype.min = function() 
{
	var min = this[0];
	var len = this.length;
	for (var i = 1; i < len; i++) 
	{
		if (this[i] < min)
		{
			min = this[i];
		}
	}
	return min;
}

function jmol_array_from_array(data)
{
	var result = "{";
	for(var i=0; i<data.length; i++)
	{
		if(i==0)
		{
			result = result + data[i][1];
			//window.console.log(data[i]);
		}
		else
		{
			result = result + "," + data[i][1];
		}
	}
	
	result = result + "}";
	
	return result;
}

//Parse a reduced mutation code ($Chain$ResidueIndex$ResidueCode$MutationCode)
//into a four element array.
function parseMutationCode(code)
{
	//Some pdbs may not have a chain id, so we have to check if one is present
	chain = code.substring(0,1);
	if((parseInt(chain) == 0) || isNaN(parseInt(chain)))
	{
		residueNumber = code.substring(1, code.length - 2);
		residueCode = code.substring(code.length -2, code.length - 1);
		mutationCode = code.substring(code.length -1, code.length);
	}
	else
	{
		chain = "";
		residueNumber = code.substring(0, code.length - 2);
		residueCode = code.substring(code.length -2, code.length - 1);
		mutationCode = code.substring(code.length -1, code.length);
	}
	
	return new Array(chain, residueNumber, residueCode,  mutationCode);
}

//Function connected to the color structure radio button array
//value is a string indicating how to colour the structure
//the valid values are structure, Scan, Binding, Stability.
//Also sets the min and max scale values.
function colorStructure(value)
{
	if(value == "structure")
	{
		jmolScript('select {*}');
		jmolScript('color ribbon structure;');
		HideContent('Scale');
	}
	else
	{	
		var values = new Array();
		
		//jmolScript('{*}.temperature=0.01');
		jmolScript('color ribbon None');

		list = '(';
		var data = dataForCalculation(value);
		for(var i=0; i<data.length; i++)
		{
			values[i] = parseFloat(data[i][1]);
			//Amazingly you can't set temperature to 0 in Jmol
			//Have to use some small non-zero number
			if(values[i] == 0.0)
			{
				data[i][1] = '0.0000001';
			}
			var script = '.CA}.temperature=';
			var codeData = parseMutationCode(data[i][0]);
			selectionString =  codeData[1];
			if(codeData[0] != "")
	   		{
				selectionString = selectionString + ':' + codeData[0];
	   		}
			script = '{'+ selectionString + script + data[i][1];
			jmolScript(script);
			list = list + selectionString + '.CA';
			if(i != data.length - 1)
				list = list + ',';
		}
		list = list + ')';
		command = 'select ' + list + ';set rangeSelected on; color ribbon relativeTemperature';
		jmolScript(command);
		
		//Set the scale values
		var min = values.min();
		var max = values.max();
		var mid = (max + min)*0.5;
		//Round to three decimal places
		mid = Math.round(mid*Math.pow(10,3))/Math.pow(10,3);
		document.getElementById('maxvalue').innerHTML = max;
		document.getElementById('minvalue').innerHTML = min;
		document.getElementById('midvalue').innerHTML = mid;
		ShowContent('Scale');
	}	
}

//Taken from http://www.codeproject.com/KB/miscctrl/JS_Inspect_Object.aspx
function inspect(obj, maxLevels, level)
{
	var str = '', type, msg;
	
	// Start Input Validations
	// Don't touch, we start iterating at level zero
	if(level == null)  level = 0;
	
	// At least you want to show the first level
	if(maxLevels == null) maxLevels = 1;
	if(maxLevels < 1)     
		return '<font color="red">Error: Levels number must be > 0</font>';
	
	// We start with a non null object
	if(obj == null)
		return '<font color="red">Error: Object <b>NULL</b></font>';
	// End Input Validations
	
	// Each Iteration must be indented
	str += '<ul>';
	
	// Start iterations for all objects in obj
	for(property in obj)
	{
		try
		{
			// Show "property" and "type property"
			type =  typeof(obj[property]);
			str += '<li>(' + type + ') ' + property + 
			( (obj[property]==null)?(': <b>null</b>'):('')) + '</li>';
			
			// We keep iterating if this property is an Object, non null
			// and we are inside the required number of levels
			if((type == 'object') && (obj[property] != null) && (level+1 < maxLevels))
				str += inspect(obj[property], maxLevels, level+1);
		}
		catch(err)
		{
			// Is there some properties in obj we can't access? Print it red.
			if(typeof(err) == 'string') msg = err;
			else if(err.message)        msg = err.message;
			else if(err.description)    msg = err.description;
			else                        msg = 'Unknown';
			
			str += '<li><font color="red">(Error) ' + property + ': ' + msg +'</font></li>';
		}
	}
	
	// Close indent
	str += '</ul>';
	
	return str;
}

//--></script>

