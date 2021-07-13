import argparse
from dataclasses import dataclass
import json
import logging
import multiprocessing
import os
import subprocess
import sys
from typing import List, Dict, Tuple

import networkx as nx


@dataclass(frozen=True, eq=True)
class Compile(object):
    """
    Data class that wraps the json Ninja output.
    """
    directory: str
    command: str
    file: str

    @property
    def include_command(self) -> List[str]:
        """
        Returns the gcc command that will output the include hierarchy.
        :return: List of strings ready to be used in a subprocess.call.
        """
        tmp = self.command.split()
        tmp.insert(1, "-E")  # gcc postprocessing option to stop after preprocessing.
        tmp.insert(1, "-H")  # gcc postprocessing option to get the included headers
        return tmp


def get_compile_commands(ninja_compile_command: str) -> Tuple[Compile]:
    """
    Converts Ninja's json into a list of Compile instances.
    :param ninja_compile_command: absolute path the the ninja `compile_commands.json` file
    :return: The tuple of compilation.
    """
    with open(ninja_compile_command, "r") as f:
        raw = json.load(f)
    return tuple(map(
        lambda j: Compile(j["directory"], j["command"], j["file"]),
        raw
    ))


def recursive_populate(edges: Tuple[str, str], root_indent: int, root_path: str, sub_graph_info: Tuple[int, str]) -> None:
    """
    Constructs the edges and populate the list.
    :param edges: List the will contain the pairs of string.
    :param root_indent: The gcc log indents more deep graph nodes.
    :param root_path: One node/file/root includes other nodes/files/children.
    :param sub_graph_info: List of (child indent, child path), respecting the gcc log.
    :return: None
    """
    for index, (indent, path) in enumerate(sub_graph_info):
        if indent == (root_indent + 1):
            edges.append((root_path, path))
            recursive_populate(edges, indent, path, sub_graph_info[index + 1:])
        if indent <= root_indent:  # This is not our graph branch anymore
            return


def find_cycles(outputs: Dict[Compile, str], exclude_dirs: List[str]) -> int:
    """
    Builds the graph of include dependencies and search for cycles.
    Cycles are printed on screen.
    :param exclude_dirs: List of directories we do not want to consider for the cycle analysis.
    :param outputs: For each Compile instance, the gcc log with the include graph (to be parsed).
    :return: The number of cycles
    """
    edges = []
    for cc, stderr in outputs.items():
        # For each compilation log, contains the depth and name for each included file.
        sub_graph_info = []
        logging.info("Populating include graph for %s." % cc.file)
        for line in filter(None, stderr.split('\n')):
            if line.startswith('.'):
                indent, path = line.split()
                # Making all the paths absolute
                if not os.path.isabs(path):
                    path = os.path.normpath(os.path.join(cc.directory, path))
                sub_graph_info.append((len(indent), path))
            else:
                break  # graph lines are on top, we stop at first non graph line.
        recursive_populate(edges, 0, cc.file, sub_graph_info)

    # Using dedicated graph lib to find cycles.
    def edge_predicate(edge: Tuple[str, str]) -> bool:
        for node in edge:
            for e in exclude_dirs:
                if node.startswith(e):
                    return False
        return True

    graph = nx.DiGraph()
    graph.add_edges_from(filter(edge_predicate, edges))
    cycles = tuple(nx.simple_cycles(graph))
    logging.info("%s cycle(s) found" % len(cycles))
    for cycle in cycles:
        print(cycle)
    return len(cycles)


def get_gcc_output(cc: Compile) -> str:
    """
    Gets the log that will be parsed into dependency graph.
    :param cc: The Compile command instance
    :return: The stderr

    sys.exit(-1) in case of error.
    """
    try:
        logging.info("Extracting include information for " + cc.file)
        process = subprocess.run(cc.include_command, capture_output=True, cwd=cc.directory, check=True)
        return process.stderr.decode('utf-8')  # output is on stderr.
    except subprocess.CalledProcessError as e:
        logging.error("Could not run " + " ".join(cc.include_command), exc_info=e)
        sys.exit(-1)


def get_gcc_outputs(ninja_compile_command: str) -> Dict[Compile, str]:
    """
    Given the list of compile commands from Ninja, gets the logs containing the included files information.
    :param ninja_compile_command: absolute path the the ninja `compile_commands.json` file
    :return: A Compile to log dictionary.
    """
    compile_commands = get_compile_commands(ninja_compile_command)
    # multiprocessing.cpu_count does not consider some limitations, si I use os.sched_getaffinity
    with multiprocessing.Pool(len(os.sched_getaffinity(0))) as pool:
        # get_gcc_output write to stdout in parallel using the logging module
        # I did not see any any interweaving so...
        logs = pool.map(get_gcc_output, compile_commands)
    return {cc: log for cc, log in zip(compile_commands, logs)}


def parse(cli_args: List[str]):
    """
    Parse the command line arguments and return the corresponding structure.
    :param cli_args: The list of arguments (as strings)
    :return: The struct
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compile-command', help="Ninja's compile_command.json file.",
                        type=str, dest="ninja_compile_command")
    parser.add_argument('-e', '--exclude-dir', default=[], action='append', dest="exclude_dirs",
                        help="Exclude a directory for cycle analysis. Can be used multiple times.")
    args = parser.parse_args(cli_args)
    return args


def main():
    """
    Exits on -1 in case of internal error. Otherwise exits on number of cycles found.
    No cycle found will therefore exit on 0.
    :return: None
    """
    logging.basicConfig(format='[%(asctime)s][%(levelname)s] %(message)s', level=logging.INFO)
    args = parse(sys.argv[1:])
    logging.info("Starting cycle detection process.")
    outputs = get_gcc_outputs(args.ninja_compile_command)
    sys.exit(find_cycles(outputs, args.exclude_dirs))


if __name__ == "__main__":
    main()
