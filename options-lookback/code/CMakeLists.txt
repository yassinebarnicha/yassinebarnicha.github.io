cmake_minimum_required(VERSION 3.16)
project(UnifiedLookback LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ==============================
# Create ONE shared DLL
# ==============================

add_library(UnifiedLookback SHARED

    # ==========================
    # BGK SOURCES
    # ==========================
    Chivet/src/AnalyticLookback.cpp
    Chivet/src/MarketParameters.cpp
    Chivet/src/MonteCarloSimulator.cpp
    Chivet/src/OptionParameters.cpp
    Chivet/src/Pricer.cpp
    Chivet/src/LookbackDLL.cpp

    # ==========================
    # (Conditional Law of max)
    # ==========================
    Barnicha/src/analytic_engine.cpp
    Barnicha/src/analytic_lookback.cpp
    Barnicha/src/bump_and_revalue_greeks.cpp
    Barnicha/src/lookback_option.cpp
    Barnicha/src/lookback_payoff_evaluator.cpp
    Barnicha/src/math.cpp
    Barnicha/src/monte_carlo_engine.cpp
    Barnicha/src/one_step_bridge_sampler.cpp
    Barnicha/src/pricing_engine.cpp
    Barnicha/src/random.cpp
    Barnicha/src/lookback_dll.cpp

    # ==========================
    # Unified export layer
    # ==========================
    unified/unified_dll.cpp
)

# ==============================
# Include directories
# ==============================

target_include_directories(UnifiedLookback PRIVATE
    Chivet/include
    Barnicha/include
)

# ==============================
# Remove lib prefix on Windows
# ==============================

set_target_properties(UnifiedLookback PROPERTIES
    PREFIX ""
)
