<?php
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

	//Contains functions for formatting generic interface elements
	
	//Sets the attributes of a datagrid.
	function set_datagrid_attributes($renderer)
	{
		$renderer->setTableAttribute("class", "results");
		$renderer->setTableAttribute("border", 1);
		$renderer->setTableHeaderAttributes(array("class"=>"results"));
		$renderer->setTableHeaderAttributes(array("style"=>"'text-wrap: nowrap'"));
		$renderer->setTableEvenRowAttributes(array("class"=>"results"));
		$renderer->setTableOddRowAttributes(array("class"=>"results"));
	}
?>
