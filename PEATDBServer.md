# PEAT Server Setup #

PEATDB allows for shared, multi-user network access to a remote database held on a server. We can utilise any storage solution provided for in ZODB. The default ZODB server-side storage is called ClientStorage and is provided by ZEO. ZEO (Zope Enterprise Objects), extends the ZODB to provide access to objects over a network.
However the default remote storage used by PEAT is RelStorage, which uses a relational database instead. This allows us to use the better security features provided by an RDBMS like MySQL, which is what we use. Supposedly RelStorage also handles high concurrency better than ZEO, though this generally won't be an issue in PEAT unless a very large number of users are reading and writing to the database at once.
Note: RelStorage simply uses the RDBMS machinery to store the object database, it is not in itself a normal relational database - i.e. we can't use mysql to search our DB.

## Setting up a relstorage server (MySQL backend) ##

RelStorage is a python library that allows one to use MySQL to store the object database, giving the security and reliablity of MySQL. Therefore this is essentially no different to setting up a MySQL server. We recommend phpMyAdmin for web based admin of the server, making it easy to add users etc. There are instructions for doing this on each operating system and rather than provide them here, we recommend a web search for that. When MySQL is setup you just add a database for each peat project and add user/passwords and give them access to the databases you want.
Finally you will need to allow outside access to the port that MySQL is running on if you want users to be able to remotely access your database. You can change the MySQL port in the configuration file too if needed. The default port is 3306. We use port 8080 on our own server. Note that because the database is open to the outside world there are security issues and you should carefully consider what data you put there. Access can be granted per IP address to individual DBs. See http://www.cyberciti.biz/tips/how-do-i-enable-remote-access-to-mysql-database-server.html

### relstorage migration to 1.5 ###

If your database was created using relstorage 1.4 and the relstorage client version has changed to 1.5 you will need to run the following MySQL code to migrate the database. This will not change any of your data:

```
    ALTER TABLE object_state ADD COLUMN state_size BIGINT AFTER md5;
    UPDATE object_state SET state_size = COALESCE(LENGTH(state), 0);
    ALTER TABLE object_state MODIFY state_size BIGINT NOT NULL AFTER md5;
```

## Server Configuration ##

You may need to add the following lines in my.cnf (the mysql configuration file):
```
[mysqld]

port = 8080 #can be any valid port
skip-external-locking
skip-name-resolve #avoids certain connection problems from remote clients
bind-address	= your ip address #usually required
```

### Granting read-only access ###
This can be done by only allowing users SELECT privileges on the database in question. In that case commits will fail if they try to save to the DB.

### Encrypted Traffic ###
Is supported in MySQL by using ssl connections. We have not tested this comprehensively, but it will be supported in PEAT.

## Setting up a ZEO Server ##

We do not recommend this method as it doesn't have any particular advantage over using the default method described above. Also, the documentation is scarce on ZEO.
Setup is as follows:
Install ZODB if you have not already
Make a folder on the server where you wish to keep the database files
Create a configuration file. We provide a default configuration file for ZEO, that is included in the peat distribution. It can also be viewed here. You will need to change the directory name and can add multiple storages.
Use zeoctl to run the server. You can test the configuration file by typing zeoctl -C zeo.conf, this will tell you if your setup is ok.
For linux systems, an init startup script is available here: Zeo.init. You can place this in the appropriate place on your system (for example on fedora this is /etc/init.d
Once you have started the server you can test by trying to connect to it from peat.
ZEO runs on port 8090 by default but this can be changed.
NB: Users will need to choose the appropriate backend option when connecting, so that the application knows the server is ZEO.

## Request Hosting of a project ##

If you don't want the trouble of setting up your own server and want to try PEATDB, send us an email and we will gladly host the project for you.

## Links ##

http://pypi.python.org/pypi/RelStorage/1.4.0b3

http://docs.zope.org/zodb/documentation/guide/zeo.html