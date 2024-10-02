# This file is part of BenchExec, a framework for reliable benchmarking:
# https://github.com/sosy-lab/benchexec
#
# SPDX-FileCopyrightText: 2007-2020 Dirk Beyer <https://www.sosy-lab.org>
#
# SPDX-License-Identifier: Apache-2.0

import logging
from xml.etree import ElementTree

from benchexec.tools.sv_benchmarks_util import get_data_model_from_task, ILP32, LP64
import benchexec.tools.template
import benchexec.result as result
import re
from pathlib import Path

import json 


class Tool(benchexec.tools.template.BaseTool2):
    """
    Internal tool info module for Modest
    """

    def executable(self, tool_locator):
        return tool_locator.find_executable("modest")

    def version(self, executable):
        return self._version_from_tool(executable)

    def name(self):
        return "Modest"

    def cmdline(self, executable, options, task, rlimits):
        model_file = task.input_files_or_identifier[0]
        print(model_file)
        print(options)
        options.insert(3, model_file)
        # No property file in modest
        property_tag = task.options.get("property_tag")
        if task.options:
            options = options + ["--props", property_tag]
        
        constants = task.options.get("constants")
        if constants:
            options += ["-E", constants]


        self.options = options

        return [executable] + options


    def analyzeValues(self, values:list[float]):
        avg = sum(values) / len(values)
        var = sum(map(lambda x: (x - avg)**2, values)) / len(values)
        maxv = max(values)
        minv = min(values)
        if avg == 0:
            cv = 0
        else:
            cv = var / avg
        return avg, var, maxv, minv, cv
    
    
    def get_value_from_output(self, output, identifier):
        if identifier == "probability":        
            value_regex = re.compile(r"robability: (.*)")
            for line in reversed(output):
                line = line.strip()
                value_match = value_regex.search(line)
                if value_match:                
                    return value_match.group(1)
        elif "distribution" in identifier:
             for line in reversed(output):
                 if "Results per scheduler" in line:
                    data = {}
                    values = list(map(float, re.findall(r'\(\d+\,\s([\d\.]+)\)', line)))
                    data["values"] = values
                    avg, var, maxv, minv, cv = self.analyzeValues(values)
                    if identifier == "distribution-mean":
                        return avg
                    elif identifier == "distribution-var":
                        return var
                    elif identifier == "distribution-maxv":
                        return maxv
                    elif identifier == "distribution-minv":
                        return minv
                    elif identifier == "distribution-coeff-var":
                        return cv
                    else:                    
                        data["mean"] = avg
                        data["variance"] = var
                        data["maxv"] = maxv
                        data["minv"] = minv
                        data["cv"] = cv
                        return json.dumps(data)
        return "NA"


