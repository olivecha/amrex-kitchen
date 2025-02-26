import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="amr_kitchen",
    version="1.0.0",
    author="Olivier Chabot",
    description=("A Toolbox to manipulate AMReX plotfiles"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/olivecha/mandoline",
    project_urls={
        "Bug Tracker": "https://github.com/olivecha/mandoline/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["requests"],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "chef = amr_kitchen.chef.cli:main",
            "chk2plt = amr_kitchen.chk2plt.cli:main",
            "colander = amr_kitchen.colander.cli:main",
            "combine = amr_kitchen.combine.cli:main",
            "mandoline = amr_kitchen.mandoline.cli:main",
            "menu = amr_kitchen.menu.cli:main",
            "pestle = amr_kitchen.pestle.cli:main",
            "taste = amr_kitchen.taste.cli:main",
            "whip = amr_kitchen.whip.cli:main",
            "minuterie = amr_kitchen.minuterie:main",
            "marinate = amr_kitchen.marinate:main",
        ]
    }
)
