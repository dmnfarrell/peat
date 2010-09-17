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


//Global variable referencing the javascript object that
//slides in and out the graphs
var collapsible;
var graph;
var data;

// Load the Visualization API and the piechart package.
google.load('visualization', '1', {'packages':['piechart','columnchart']});
// Set a callback to run when the Google Visualization API is loaded.
google.setOnLoadCallback(CreateGraph);


function HideContent(d) 
{
	document.getElementById(d).style.display = "none";
}

function ShowContent(d) 
{
	document.getElementById(d).style.display = "block";
}

function ReverseDisplay(d) 
{
	if(document.getElementById(d).style.display == "none") 
	{ 
		document.getElementById(d).style.display = "block"; 
	}
	else
	{ 
		document.getElementById(d).style.display = "none"; 
	}
}

//Sets up the variables necessary for showing and hiding the results graphs.
//Initialises the global variable 'collapsible' which is the javascript object
//that animates the opening and closing of the graph image.
//Also hides the image initially
function InitialiseGraph()
{
	collapsible = new Fx.Slide($('resultsGraph'), {
				   duration: 500,
				   transition: Fx.Transitions.linear
				   });
	HideGraphBox();			   
}

function CreateGraph()
{
	graph = new google.visualization.ColumnChart(document.getElementById('resultsGraph'));
}

//This function is called when one graph is displayed and another is requested
function HandleGraphChange()
{
	graph.draw(data, {height: 400, legend: 'none', fontSize: 12, title: document.getElementById('resultsGraph').name});
	//Slide the image out and remove this handler from the collapsible object
	collapsible.slideIn();
	collapsible.removeEvent('complete', HandleGraphChange);
}

function SetGraphData(resultType) 
{
	// Create our data table.
	var calculationData = dataForCalculation(resultType);
	data = new google.visualization.DataTable();
	
	//Stability and Binding calcs just report Mutation and Total
	//For pKa calls each mutation will posisble have results for a number of ionisable residues
	if(resultType == "Stability" || resultType == "Binding")
	{
		data.addColumn('string', 'Mutation');
		data.addColumn('number', 'Free-Energy Difference');

		for(var i=0; i < calculationData.length; i++)
		{
			calculationData[i][1] = parseFloat(calculationData[i][1]);
			data.addRow(calculationData[i]);
		}
	}
	else if(resultType == "Scan")
	{
		data.addColumn('string', 'Mutation');
		for(i=1; i<calculationData[0].length; i++)
		{
			data.addColumn('number', calculationData[0][i]);
		}

		for(var i=1; i < calculationData.length; i++)
		{
			for(j=1; j<calculationData[i].length; j++)
			{
				calculationData[i][j] = parseFloat(calculationData[i][j]);
			}
			data.addRow(calculationData[i]);
		}
	}
	else if(resultType == "Modelling")
	{
		data.addColumn('string', 'Mutation');
		data.addColumn('number', 'Bump Score');
		
		for(var i=0; i < calculationData.length; i++)
		{
			calculationData[i][1] = parseFloat(calculationData[i][1]);
			data.addRow(calculationData[i]);
		}
	}	
}

//This function shows and hides the results graphs on the result page
//There are three possibilities
//1 - Open a graph when none is displayed
//2 - Close the displayed graph
//3 - Open a new graph when one is already displayed
//resultType is one of Stability, Binding or deltapKa
function ShowResultsGraph(resultType)
{
	//Check is a new graph requested i.e. it the name different to the resultsGraph name attribute	
	if(document.getElementById('resultsGraph').name != resultType)
	{
		SetGraphData(resultType);
		document.getElementById('resultsGraph').name = resultType;
		//If an graph is already displayed we close it, replace it, and open it.
		//Otherwise we just replace the current chart and display it
		if(collapsible.open == 1)
		{
			//Set the HandleGraphChange function to be called after
			//we slide the graph out
			collapsible.addEvent('complete', HandleGraphChange);
			//Set the new link for the graph - This will also be used
			//by HandleGraphChange to obtain the new url
			//document.getElementById('resultsGraphLink').href = url;
			//Slide the graph out - when this is done HandleGraphChange will execute.
			collapsible.slideOut();
		}
		else
		{
			graph.draw(data, {height: 400, legend: 'none', fontSize: 12, title: resultType, titleY: 'kj/mol'});
			//document.getElementById('resultsGraphLink').href = url;
		        collapsible.toggle();			       
		}
	}
	else
	{
		//Its the same graph as was last selected - just toggle
		collapsible.toggle();
	}
}

function HideGraphBox()
{	
	collapsible.hide();			 
	ShowContent('resultsGraph');
}

function HideElement(id)
{
	element = new Fx.Slide($(id), {
				   duration: 500,
				   transition: Fx.Transitions.linear
				   });  
	element.hide()				
}

function ToggleElement(id)
{	
	element = new Fx.Slide($(id), {
				   duration: 250,
				   transition: Fx.Transitions.linear
				   });   
	element.toggle();		 
}


//--></script>

