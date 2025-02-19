# ECMWF / CDS
Follow instructions [here](https://cds.climate.copernicus.eu/how-to-api) to setup the CDS API personal access token.

Once logged in, copy the API Token displayed in your profile below to the file $HOME/.cdsapirc (in your Unix/Linux environment)

```
url: https://cds.climate.copernicus.eu/api
key: <PERSONAL-ACCESS-TOKEN>
```

In your environment/virtual environment download the packages:
```
pip install "cdsapi>=0.7.4" eccodes pygrib matplotlib
```

##

