import json
import math
import subprocess
import sys

# ex run: python rodeo/benchmark.py doublepir rodeo/workload.json rodeo/doublepir.json


def run_benchmark(scheme: str, num_items: int, item_size_bits: int, num_clients: int):
    log_n = math.ceil(math.log(num_items, 2))
    d = item_size_bits

    bench_name = "SimplePirSingle"
    if scheme == "doublepir":
        bench_name = "DoublePirSingle"

    out_json_filename = "pir/report.json"

    cmd = f"LOG_N={log_n} D={d} go test -bench {bench_name} -timeout 0 -benchtime=1x -run=^$"
    if num_clients > 1:
        cmd = "DO_K_TEST=1 K_CCB=" + str(num_clients) + " " + cmd

    # run commmand and get output
    subprocess.run(cmd, shell=True, cwd="./pir")

    # read output json file
    json_str = None
    with open(out_json_filename, "r") as fh:
        json_str = fh.read()

    # print json
    data = json.loads(json_str)

    # delete json file
    subprocess.check_output(f"rm {out_json_filename}", shell=True)

    return data


def run_benchmarks(scheme: str, workload_file: str, output_json_file: str):
    workload_json = None
    with open(workload_file, "r") as fh:
        workload_json = json.load(fh)

    template_json_file = "rodeo/rodeo.json"
    merged_json = None
    with open(template_json_file, "r") as fh:
        merged_json = json.load(fh)

    results = []
    for scenario in workload_json["workloads"]:
        num_items = int(scenario["db"]["numItems"])
        item_size_bits = int(scenario["db"]["itemSizeBits"])
        num_clients = (
            int(scenario["clients"]["numClients"]) if "clients" in scenario else 1
        )
        measurement = run_benchmark(scheme, num_items, item_size_bits, num_clients)
        result = {"scenario": scenario, "measurement": measurement}
        results.append(result)

        # merge and write
        merged_json["results"] = results
        merged_json["scheme"]["variant"] = scheme
        with open(output_json_file, "w") as fh:
            json.dump(merged_json, fh, indent=2)


if __name__ == "__main__":
    scheme = sys.argv[1]
    workload_file = sys.argv[2]
    output_json_file = sys.argv[3]
    run_benchmarks(scheme, workload_file, output_json_file)
