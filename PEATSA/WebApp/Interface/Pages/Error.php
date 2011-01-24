<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 TRANSITIONAL//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>

	<head>
		<title>PEAT-SA (Beta) - Job Error</title>
	<link href="SiteWideStyle.css" rel="stylesheet" type="text/css">
	</head>
	<body>
	<div id="container">
	<h1><span>Houston, we have a problem ...</span></h1>
	<h3>An error was detected</h3>
	
	<?php 
		include 'PageComponents.php';
		add_navigation_pane();		

		//If some data was supplied go through it and see what it is
		echo '<div id="content">';
		if(count($_GET) != 0)
		{	
			echo '<p>';
			if(array_key_exists("jobId", $_GET))
				echo "Job Id: ".$_GET['jobId']."<br>";
			
			if(array_key_exists("domain", $_GET))
				echo "The error occurred in the ".$_GET['domain']." - Here's the error info:<br><br>";
		
			if(array_key_exists("description", $_GET))
				echo $_GET['description']."<br>";
			else
				echo "Unknown Error.<br>";
	
			if(array_key_exists("detailedDescription", $_GET))
				echo stripslashes($_GET['detailedDescription'])."<br><br>";	
		
			echo '</p>';

			if(array_key_exists("recoverySuggestion", $_GET))
				echo "<h3>What should I do?</h3>";
				echo "<p>";
				echo stripslashes($_GET['recoverySuggestion']);
				echo '</p>';

			echo '<p>If the above information makes no sense to you please contact us. ';
			echo "We'll be happy to help you out.<br>Also check the FAQ - there may be more information on your error there.</p>";
		}
		else
		{
			echo "<p>Strangely no error data was sent to this page.<br>";
			echo "Possibly you navigated directly here?</p>";
		}

		echo "</div>";

		add_colophon();
	?>
	</div>
	</body>
</html>
