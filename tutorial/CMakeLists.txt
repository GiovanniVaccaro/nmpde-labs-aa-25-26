cmake_minimum_required(VERSION 3.12.0)
project(tutorial-lab LANGUAGES CXX C)

# Carica la configurazione comune del corso (fondamentale per usare deal_ii_setup_target).
# Questo assume che la cartella "common" sia nella directory superiore, come nel lab02.
include(../common/cmake-common.cmake)

# --- Tutorial 01 ---
# Modifica "tutorial-01.cpp" se il tuo file ha un nome leggermente diverso (es. tutorial01.cpp)
add_executable(tutorial-01 src/tutorial-01.cpp)
deal_ii_setup_target(tutorial-01)

# --- Tutorial 02 ---
add_executable(tutorial-02 src/tutorial-02.cpp)
deal_ii_setup_target(tutorial-02)

# --- Tutorial 03 ---
add_executable(tutorial-03 src/tutorial-03.cpp)
deal_ii_setup_target(tutorial-03)

# --- Tutorial 04 ---
add_executable(tutorial-04 src/tutorial-04.cpp)
deal_ii_setup_target(tutorial-04)

# --- Tutorial 05 ---
add_executable(tutorial-05 src/tutorial-05.cpp)
deal_ii_setup_target(tutorial-05)