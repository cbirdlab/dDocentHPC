/*
 * This mpi_parallel implemented 0.1% functionality of the mighty GNU parallel. 
 *
 * It starts process with MPI instead of on local node or through ssh
 * checkout GNU parallel for more power: 
 *   https://www.gnu.org/software/parallel/
 * Please cite GNU parallel to support Open Source community:
 *   http://git.savannah.gnu.org/cgit/parallel.git/tree/doc/citation-notice-faq.txt
 *
 * Usage: mpi_parallel command [a] {} [b]....
 * read stdin line by line and replace {} with each line
 *
 * environment variables for settings
 * MPI_PARALLEL_VERBOSE        - display debug information 
 * MPI_PARALLEL_IGNORE_ERROR   - if any command exit without return 0, MPI do not terminal entirely
 *
 * environment variables set for child process
 * MPI_PARALLEL_HOST           - hostname
 * MPI_PARALLEL_RANK           - MPI rank 
 * MPI_PARALLEL_INJECT_*       - all environment variables starts with MPI_PARALLEL_INJECT 
 *                               will be injected right before command launching
 *                               format: MPI_PARALLEL_INJECT_foo=bar
 *                                       foo=bar will be injected
 * example:
 * 
 * $ cat input
 * test1
 * test2
 * test3
 *
 * $ export MPI_PARALLEL_INJECT_TEST_VAR="hello"
 *
 * $ cat input | mpirun -n 3 ./mpi_parallel 'echo $TEST_VAR {} from $MPI_PARALLEL_HOST $MPI_PARALLEL_RANK'
 * hello test3 from wahab-02 2
 * hello test1 from wahab-02 0
 * hello test2 from wahab-02 1
 */
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <mpi.h>
#include <stdlib.h>

#define STDIN_INPUT_LINE_TAG  0x010

#define INFO(x) do {  \
  if (MPI_PARALLEL_VERBOSE) \
  std::cout << __func__ << " " << __LINE__ << ": " << x << std::endl;  \
} while(0)

#define ERROR(x, y) do {  \
  std::cerr << x << std::endl;  \
  MPI_Abort(MPI_COMM_WORLD, y); \
  return y; \
} while(0)

bool MPI_PARALLEL_VERBOSE      = false;
bool MPI_PARALLEL_IGNORE_ERROR = false;

extern char **environ;

std::tuple<int, int, std::string> mpi_kickoff(){
  int size;
  int rank;

  char hostname[MPI_MAX_PROCESSOR_NAME];
  int  hostname_len;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(hostname, &hostname_len);

  return std::make_tuple(size, rank, std::string(hostname));
}

void command_arguments_scatter(int world_size) {
  using namespace std;
  vector<string> inputs;
  string line;

  while (getline(cin, line)) {
    inputs.push_back(line);
  }

  if (inputs.size() > world_size) {
    stringstream ss;
    ss <<  "Insufficient MPI job size: MPI job(" << world_size << ") input size("  << inputs.size() << ")" ;
    throw ss.str();
  }

  int i = 0;

  for (auto line : inputs) {
    MPI_Send(line.c_str(), line.size() + 1, MPI_CHAR, i++, STDIN_INPUT_LINE_TAG, MPI_COMM_WORLD);
  }

  while(i < world_size) {
    MPI_Send("", 1, MPI_CHAR, i++, STDIN_INPUT_LINE_TAG, MPI_COMM_WORLD);
  }
}

std::string command_arguments_recv() {
  char input_buffer[512];
  MPI_Recv(&input_buffer, 512, MPI_CHAR, MPI_ANY_SOURCE, STDIN_INPUT_LINE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return std::string(input_buffer);
}

std::string command_arguments_replace(int argc, char **argv, std::string input) {
  using namespace std;
  stringstream ss;
  for (int i = 1; i < argc; ++i) {
    ss << argv[i] << " ";
  }
  auto command_template = ss.str();
  auto found = command_template.find("{}");
  while (found!=string::npos){
    command_template = command_template.replace(found, 2, input);
    found = command_template.find("{}");
  }

  return command_template;
}

void environment_variables_inject(std::string hostname, int world_rank) {
  using namespace std;

  string target = "MPI_PARALLEL_INJECT_";

  for (char **raw_env = environ; *raw_env; ++raw_env) {
    auto env = string(*raw_env);
    if (env.compare(0, target.size(), target) != 0)
      continue;

    auto pos = env.find("=");
    if (pos==string::npos)
      continue;
    auto env_name = env.substr(target.size(), pos-target.size());
    auto env_value = env.substr(pos + 1, string::npos);

    setenv(env_name.c_str(), env_value.c_str(), 1);
    INFO("stdenv: " << env_name << "=" << env_value);
  }

  setenv("MPI_PARALLEL_HOST", hostname.c_str(), 1);
  setenv("MPI_PARALLEL_RANK", to_string(world_rank).c_str(), 1);
}

int main(int argc, char **argv) {
  using namespace std;

  try {
    //c++ 17?
    //auto [world_size, world_rank, hostname] = MPI_kickoff();
    //
    int world_size;
    int world_rank;
    string hostname;
    tie(world_size, world_rank, hostname) = mpi_kickoff();


    if (getenv("MPI_PARALLEL_VERBOSE"))
      MPI_PARALLEL_VERBOSE = true;

    if (getenv("MPI_PARALLEL_IGNORE_ERROR"))
      MPI_PARALLEL_IGNORE_ERROR = true;

    if (world_rank == 0) {
      command_arguments_scatter(world_size);
    }

    auto arguments = command_arguments_recv();

    if (arguments.size() == 0) {
      INFO(hostname << " " << world_rank << " out of " << world_size << " empty input, exiting normally");
      MPI_Finalize();
      return 0;
    }
    auto command = command_arguments_replace(argc, argv, arguments);
    INFO(hostname << " " << world_rank << " out of " << world_size << " " << command);
    environment_variables_inject(hostname, world_rank);

    int exit_code = system(command.c_str());

    if (exit_code != 0 && !MPI_PARALLEL_IGNORE_ERROR) {
      ERROR(hostname << " " << world_rank << " error code: " << exit_code <<" when running: '" << command << "'", exit_code);
    }

    MPI_Finalize();
  } catch (string msg){
    ERROR(msg, -1);
  }

  return 0;
}
