import setuptools


with open("README.md") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EllipticGenus",
    version="1.0.0",
    author="Kenta Kobayashi, Akihito Nakamura, Kazushi Ueda",
    author_email="kenta.topos@gmail.com",
    license="GPL2+",
    description="Computation of elliptic genera of homogeneous spaces and their complete intersections",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    zip_safe=False,
)
