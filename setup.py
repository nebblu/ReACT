from setuptools import setup, find_packages
from distutils.command.install import install
from distutils.command.build import build
from distutils.command.clean import clean

import subprocess

import os
import sys
import shutil

LIBS = ["libreact_wrapper.so",]
HERE = os.path.abspath(os.path.dirname(__file__))

def compile_library(env):
    subprocess.check_call(["make"], env=env, cwd="pyreact")
    # shutil.copy(os.path.join(here, "cl_to_xi/libwigner.so"), os.path.join(here, "tpst/libwigner.so"))

def clean_library(env={}):
    subprocess.check_call(["make", "clean"], env=env, cwd="pyreact")
    # os.remove(os.path.join(here, "tpst/libwigner.so"))

class my_build(build):
    def run(self):
        env = os.environ
        compile_library(env)
        super().run()


class my_install(install):
    def __init__(self, dist):
        install.__init__(self, dist)
        self.build_args = {}
        if self.record is None:
            self.record = "install-record.txt"

    def run(self):
        super().run()

class my_clean(clean):
    def run(self):
        clean_library()
        super().run()

setup(name = "pyreact",
      description       = "Python interface for the ReACT library",
      author            = "Tilman Troester",
      author_email      = "tilmantroester@gmail.com",
      packages = ["pyreact"],
      package_data = {"" : LIBS,},
      install_requires = ['numpy'],
      cmdclass={"install"   : my_install,
                "build"     : my_build,
                "build_ext" : my_build,
                "clean"     : my_clean},
        )

