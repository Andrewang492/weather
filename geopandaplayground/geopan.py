import geopandas
from geodatasets import get_path
import matplotlib.pyplot as plt
path_to_data = get_path("nybb")
gdf = geopandas.read_file(path_to_data)

gdf

gdf = gdf.set_index("BoroName")
gdf["area"] = gdf.area
gdf["area"]

gdf["boundary"] = gdf.boundary
gdf["centroid"] = gdf.centroid
first_point = gdf["centroid"].iloc[0]
gdf["distance"] = gdf["centroid"].distance(first_point)


gdf.plot("area", legend=True)
plt.savefig("area.png")
gdf.explore("area", legend=False)
gdf = gdf.set_geometry("centroid")
gdf.plot("area", legend=True)

ax = gdf["geometry"].plot()
gdf["centroid"].plot(ax=ax, color="black")
plt.savefig("centroid.png")

gdf = gdf.set_geometry("geometry") # reset the geometry column

#If we are interested in the convex hull of our polygons, we can access GeoDataFrame.convex_hull.

gdf["convex_hull"] = gdf.convex_hull
# saving the first plot as an axis and setting alpha (transparency) to 0.5
ax = gdf["convex_hull"].plot(alpha=0.5)
# passing the first plot and setting linewidth to 0.5
gdf["boundary"].plot(ax=ax, color="white", linewidth=0.5)
#lt.savefig("convex_hull.png")

# buffering the active geometry by 10 000 feet (geometry is already in feet)
gdf["buffered"] = gdf.buffer(10000)

# buffering the centroid geometry by 10 000 feet (geometry is already in feet)
gdf["buffered_centroid"] = gdf["centroid"].buffer(10000)
# saving the first plot as an axis and setting alpha (transparency) to 0.5
ax = gdf["buffered"].plot(alpha=0.5)
# passing the first plot as an axis to the second
gdf["buffered_centroid"].plot(ax=ax, color="red", alpha=0.5)
#plt.savefig("bufferedcentroid.png")
# passing the first plot and setting linewidth to 0.5
gdf["boundary"].plot(ax=ax, color="white", linewidth=0.5)
#plt.savefig("bufferedboundary.png")

# Geometry relations
# get the polygon of brooklyn
brooklyn = gdf.loc["Brooklyn", "geometry"]
print(type(brooklyn))
print(gdf["buffered"].intersects(brooklyn))

# check intersections by buffered centroids. Check what is entirely within.
gdf["within"] = gdf["buffered_centroid"].within(gdf)
print(gdf["within"])
gdf = gdf.set_geometry("buffered_centroid")
# using categorical plot and setting the position of the legend
ax = gdf.plot(
    "within", legend=True, categorical=True, legend_kwds={"loc": "upper left"}
)
# passing the first plot and setting linewidth to 0.5
gdf["boundary"].plot(ax=ax, color="black", linewidth=0.5)
plt.savefig("within.png")

## Projections on Coordinate Reference Systems
print(gdf.crs)
gdf = gdf.set_geometry("geometry")
boroughs_4326 = gdf.to_crs("EPSG:4326")
boroughs_4326.plot()
plt.savefig("boroughs_4326.png")