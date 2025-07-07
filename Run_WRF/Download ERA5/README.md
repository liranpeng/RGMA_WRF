1. Setup the CDS API personal access token
Here is how to setup the CDS API personal access token:
	1. If you do not have an account yet, please register.
	2. If you are not logged in, please login.
	3. Once logged in, copy the code displayed below to the file $HOME/.cdsapirc
(in your Unix/Linux environment)
Please do not copy below:
url: https://cds.climate.copernicus.eu/api
key: b0b31706-37bd-4455-bb03-68b44476a5ba

2. Install the CDS API client
The CDS API client is a Python based library. It provides support for Python 3.
You can install the CDS API client via the package management system pip, by running on Unix/Linux the command below.
$ pip install "cdsapi>=0.7.4"
Version 0.7.2 or higher is required in order to be able to use the new data stores. Use of latest version available is highly recommended.

