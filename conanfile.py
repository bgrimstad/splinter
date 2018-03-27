from   conans       import ConanFile, CMake, tools
from   conans.tools import download, unzip
import os

class SplinterConan(ConanFile):
    name            = "splinter"
    version         = "3.0"     
    description     = "Conan package for bgrimstad splinter."           
    url             = "https://github.com/bgrimstad/splinter"
    license         = "MPL"                                         
    settings        = "arch", "build_type", "compiler", "os"
    generators      = "cmake"
    options         = {"shared": [True, False]} 
    default_options = "shared=False"

    def source(self):
        zip_name = "v%s.zip" % self.version
        download ("%s/archive/%s" % (self.url, zip_name), zip_name, verify=False)
        unzip    (zip_name)
        os.unlink(zip_name)

    def build(self):
        cmake = CMake(self)
        self.run("cmake %s-%s %s" % (self.name, self.version, cmake.command_line))
        self.run("cmake --build . %s" % cmake.build_config)

    def package(self):
        include_folder = "%s-%s/include" % (self.name, self.version)       
        self.copy("*.h"  , dst="include/splinter", src=include_folder)
        self.copy("*.a"  , dst="lib", keep_path=False)
        self.copy("*.so" , dst="lib", keep_path=False)
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["splinter-3-0"] if self.options.shared else ["splinter-static-3-0"]
