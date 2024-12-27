

import sys

def read_graph(file_path):
  graph = {}
  with open(file_path, 'r') as file:
    for line in file:
      if line.startswith('#'):
        continue
      if ' ' in line:
        v1, v2 = map(int, line.strip().split())
      else:
        v1, v2 = map(int, line.strip().split(','))
      if v1 not in graph:
        graph[v1] = set()
      if v2 not in graph:
        graph[v2] = set()
      graph[v1].add(v2)
      graph[v2].add(v1)
  return graph

def write_graph(graph, file_path):
  with open(file_path, 'w') as file:
    for v in sorted(graph):
      for adj in sorted(graph[v]):
        file.write(f"{v},{adj}\n")

def main(input_file, output_file):
  """
  Converts a directed graph to an undirected graph.
  """
  graph = read_graph(input_file)
  write_graph(graph, output_file)

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python to_directed.py <input_file> <output_file>")
    sys.exit(1)
  input_file = sys.argv[1]
  output_file = sys.argv[2]
  main(input_file, output_file)