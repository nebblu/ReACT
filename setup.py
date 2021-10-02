from setuptools import setup, find_packages
from distutils.command.install import install as _install
from distutils.command.clean import clean as _clean
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop

import subprocess

import os
import sys
import shutil

LIBS = ["libreact_wrapper.so",]
HERE = os.path.abspath(os.path.dirname(__file__))

def compile_library(env):
    subprocess.check_call(["make"], env=env, cwd="pyreact")

def clean_library(env={}):
    subprocess.check_call(["make", "clean"], env=env, cwd="pyreact")

class build(_build_py):
    def run(self):
        env = os.environ
        compile_library(env)
        super().run()

class develop(_develop):
    def run(self):
        env = os.environ
        compile_library(env)
        super().run()

class install(_install):
    def __init__(self, dist):
        super().__init__(dist)
        self.build_args = {}
        if self.record is None:
            self.record = "install-record.txt"

    def run(self):
        super().run()

class clean(_clean):
    def run(self):
        clean_library()
        super().run()

setup(name = "pyreact",
      description       = "Python interface for the ReACT library",
      author            = "Tilman Troester",
      author_email      = "tilmantroester@gmail.com",
      packages = ["pyreact"],
      version='1.0.0',
      package_data = {"" : LIBS,},
      install_requires = ['numpy'],
      cmdclass={"install"   : install,
                "develop"   : develop,
                "build_py"  : build,
                "clean"     : clean},
        )
