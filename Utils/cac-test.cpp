#include "GC/ShareParty.h"
#include "GC/ShareParty.hpp"
#include "GC/ShareSecret.hpp"
#include "GC/MaliciousRepSecret.h"
#include "GC/RepPrep.h"

#include "GC/Machine.hpp"
#include "GC/Processor.hpp"
#include "GC/Program.hpp"
#include "GC/Thread.hpp"
#include "GC/ThreadMaster.hpp"
#include "GC/RepPrep.hpp"

#include "Processor/Instruction.hpp"
#include "Protocols/MaliciousRepMC.hpp"
#include "Protocols/MAC_Check_Base.hpp"
#include "Protocols/Beaver.hpp"

#include "Machines/ShamirMachine.hpp"
#include "Machines/Rep4.hpp"
#include "Machines/Rep.hpp"
#include "Protocols/ProtocolSet.h"

#include <vector>
#include "Tools/time-func.h"

using namespace std;


template<class T>
void run(char** argv);

int main(int argc, char** argv) {
    // need player number and number of players
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0]
                << "<my number: 0/1/...>"
                << endl;
        exit(1);
    }

    run<GC::MaliciousRepSecret>(argv);

    return 0;
}

template<class T>
void run(char** argv)
{
    // run 64-bit computation by default
    int n_bits = 128;

    // set up networking on localhost
    int my_number = atoi(argv[1]);
    int n_parties = 3;
    int port_base = 9999;
    Names N(my_number, n_parties, "localhost", port_base);
    CryptoPlayer P(N);

    OnlineOptions::singleton.batch_size = 20000;

    // protocol setup (domain, MAC key if needed etc)
    BinaryProtocolSetup<T> setup(P);

    // set of protocols (input, multiplication, output)
    BinaryProtocolSet<T> set(P, setup);
    auto& protocol = set.protocol;
    auto* prep = protocol.get_prep();

    Timer timer;
    timer.start();

    int n = atoi(argv[2]);
    vector<array<T, 3>> triples(n);

    for (int i = 0; i < n; i ++) {
        triples[i] = prep->get_triple(n_bits);
    }

    timer.stop();
    
    size_t total_comm = 0;
    auto comm = P.total_comm();

    for (auto each: comm) {
        total_comm += each.second.data;
    }

    cout << "Time for " << n << " triples: " << timer.elapsed() << endl;
    cout << "Comm of me is " << total_comm / (double) (1024 * 1024) << endl;
    
}