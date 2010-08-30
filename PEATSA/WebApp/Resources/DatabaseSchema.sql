--
-- Protein Engineering Analysis Tool Structure Analysis (PEATSA)
-- Copyright (C) 2010 Michael Johnston & Jens Erik Nielsen
--
-- Author: Michael Johnston
--
-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--
-- Contact information:
-- Email: Jens.Nielsen_at_gmail.com
-- Normal mail:
-- Jens Nielsen
-- SBBS, Conway Institute
-- University College Dublin
-- Dublin 4, Ireland
--
--
-- phpMyAdmin SQL Dump
-- version 3.2.3
-- http://
--
-- Host: localhost
-- Generation Time: Aug 26, 2010 at 01:52 PM
-- Server version: 5.0.88
-- PHP Version: 5.2.9

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";

--
-- Database: `PDT`
--www.phpmyadmin.net

-- --------------------------------------------------------

--
-- Table structure for table `Data`
--

CREATE TABLE IF NOT EXISTS `Data` (
  `JobID` varchar(80) NOT NULL,
  `DataSetName` varchar(80) NOT NULL default 'Output',
  `MatrixName` varchar(80) NOT NULL,
  `Content` longtext NOT NULL,
  `Size` bigint(20) unsigned NOT NULL,
  `BinData` longblob,
  PRIMARY KEY  (`JobID`,`MatrixName`),
  KEY `JobID` (`JobID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `Images`
--

CREATE TABLE IF NOT EXISTS `Images` (
  `JobID` varchar(80) NOT NULL,
  `Name` varchar(80) NOT NULL,
  `Size` bigint(20) NOT NULL,
  `Content` longblob,
  PRIMARY KEY  (`JobID`,`Name`),
  KEY `JobID` (`JobID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `Jobs`
--

CREATE TABLE IF NOT EXISTS `Jobs` (
  `JobID` varchar(80) NOT NULL,
  `State` enum('UnderConstruction','Ready','Launched','Running','Finished') NOT NULL default 'UnderConstruction',
  `Email` varchar(50) NOT NULL default 'Unknown',
  `SentMail` tinyint(1) NOT NULL default '0',
  `PDBID` varchar(30) NOT NULL,
  `Structure` longtext,
  `Ligand` longtext,
  `MutationCommand` enum('--mutation','--mutationList') default NULL,
  `MutationData` longtext,
  `Date` datetime NOT NULL,
  `Stability` enum('NotSelected','Selected','Queued','Waiting','Running','Finished') NOT NULL default 'NotSelected',
  `Scan` enum('NotSelected','Selected','Queued','Waiting','Running','Finished') NOT NULL default 'NotSelected',
  `Binding` enum('NotSelected','Selected','Queued','Waiting','Running','Finished') NOT NULL default 'NotSelected',
  `QueueStatusMessage` varchar(200) NOT NULL default 'Unknown',
  `Error` tinyint(3) unsigned NOT NULL default '0',
  `ErrorDescription` varchar(500) NOT NULL default '""',
  `DetailedDescription` varchar(1000) NOT NULL default '""',
  `Log` longtext,
  `Metadata` blob,
  PRIMARY KEY  (`JobID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `Data`
--
ALTER TABLE `Data`
  ADD CONSTRAINT `Data_ibfk_1` FOREIGN KEY (`JobID`) REFERENCES `Jobs` (`JobID`) ON DELETE CASCADE,
  ADD CONSTRAINT `Data_ibfk_2` FOREIGN KEY (`JobID`) REFERENCES `Jobs` (`JobID`) ON DELETE CASCADE;

--
-- Constraints for table `Images`
--
ALTER TABLE `Images`
  ADD CONSTRAINT `Images_ibfk_1` FOREIGN KEY (`JobID`) REFERENCES `Jobs` (`JobID`) ON DELETE CASCADE;
