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
-- Create users called 'pdt'@'localhost' and 'pdt'@'%' for the database created using DatabaseSchema.sql
-- Remeber to change 'InsertPassword' to whatever you want the password to be
--

CREATE USER 'pdt'@'localhost';
SET PASSWORD FOR 'pdt'@'localhost' = PASSWORD('InsertPassword');
GRANT USAGE ON *.* TO 'pdt'@'localhost';
GRANT INSERT, SELECT, DELETE, UPDATE ON PDT.* TO 'pdt'@'localhost';
FLUSH PRIVILEGES;

CREATE USER 'pdt'@'%';
SET PASSWORD FOR 'pdt'@'%' = PASSWORD('InsertPassword');
GRANT USAGE ON *.* TO 'pdt'@'%';
GRANT INSERT, SELECT, DELETE, UPDATE ON PDT.* TO 'pdt'@'%';
FLUSH PRIVILEGES;