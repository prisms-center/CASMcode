g++ -o casm-complete -I../../include -L/home/jgg/.local/lib -std=c++11 complete.cpp -lcasm -lboost_program_options -lboost_system -lboost_regex
g++ -o complete_test -I../../include -L/home/jgg/.local/lib -std=c++11 complete_tests.cpp -lcasm -lboost_program_options -lboost_system -lboost_regex
