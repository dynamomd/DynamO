#!/usr/bin/env python3
import os
import re
import sys
import sysconfig
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        #print("## Current working directory", os.getcwd())
        #print("Directory contents", os.listdir('.'))
        #print("Submodule contents", os.listdir('./submodule'))
        #print("## Running git submodule update", ['git', 'submodule', 'update', '--init', '--recursive'], "cwd=",ext.sourcedir)
        #subprocess.check_call(['git', 'submodule', 'update', '--init', '--recursive'], cwd=ext.sourcedir)
        #print("## Running cmake with arguments", ['cmake', ext.sourcedir] + cmake_args, "cwd=",self.build_temp, "env=",env)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        #print("## Running cmake build", ['cmake', '--build', '.'] + build_args, "cwd=",self.build_temp)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name="pydynamo",
    version="0.0.1",
    author="Marcus Bannerman",
    author_email="m.bannerman@gmail.com",
    description="PyDynamO: Controls and interfaces for the DynamO simulation environment",
    long_description=open("Readme.md", "r").read(),
    long_description_content_type="text/markdown",
    url="https://dynamomd.com",
    classifiers=[
        "Programming Language :: C++",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages('src'),
    package_dir={"":"src"},
    package_data={},
    ext_modules=[
        #CMakeExtension('pydynamo.core')
    ],
    install_requires = [
        'cmake',
        'numpy',
        'scipy',
    ],
    setup_requires = [
        "cmake",
        "ninja"
    ],
    test_suite='tests',
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
