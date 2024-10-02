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


class Tool(benchexec.tools.template.BaseTool2):
    """
    Internal tool info module for Storm
    """

    def executable(self, tool_locator):
        return tool_locator.find_executable("painful-script.sh")

    def version(self, executable):
        return self._version_from_tool(executable)

    def name(self):
        return "123Go"

    def cmdline(self, executable, options, task, rlimits):        
        model_file = task.input_files_or_identifier[0]
        options = ["--prism"] + [model_file] + options
        property_file = task.options.get("property_file")
        if property_file:
            property_file = Path(model_file).parent / property_file
        property_tag = task.options.get("property_tag")
        if task.options:
            options = options + ["--prop", property_file, property_tag]
        
        constants = task.options.get("constants")
        if constants:
            options += ["--constants", constants]

        self.options = options

        return [executable] + options

    def get_value_from_output(self, output, identifier):
        if identifier not in {"probability", "invalid_action_prediction"}:
            return "NA"
        
        if identifier == "probability":
            value_regex = re.compile(r"Result [a-zA-Z\(\)\ \-]*: (.*)")
            for line in reversed(output):
                line = line.strip()
                value_match = value_regex.search(line)
                if value_match:                
                    return value_match.group(1)
        
        if identifier == "invalid_action_prediction":
            value_regex = re.compile(r"Random-choices[a-zA-Z\(\)\ \-]*: (.*)")
            for line in reversed(output):
                line = line.strip()
                value_match = value_regex.search(line)
                if value_match:
                    return value_match.group(1)
        
        return "NO_VALUE"
