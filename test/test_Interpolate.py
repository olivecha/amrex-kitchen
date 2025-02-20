

from amr_kitchen import PlotfileCooker

pck = PlotfileCooker("../test_assets/example_plt_3d")
print(pck["temp"](0.012,0.012,0.012))
