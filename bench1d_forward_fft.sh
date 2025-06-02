#!/bin/bash
# bench1d.sh

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Setting up benchmark environment...${NC}"

# Save current settings
ORIGINAL_GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
ORIGINAL_TURBO=$(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo "N/A")
ORIGINAL_RANDOMIZE=$(cat /proc/sys/kernel/randomize_va_space)

echo "Saved original settings:"
echo "  CPU Governor: $ORIGINAL_GOVERNOR"
echo "  Turbo disabled: $ORIGINAL_TURBO"
echo "  ASLR: $ORIGINAL_RANDOMIZE"

# Function to restore settings
restore_settings() {
    echo -e "${BLUE}Restoring original settings...${NC}"
    
    # Restore CPU governor
    echo $ORIGINAL_GOVERNOR | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor > /dev/null
    echo "  CPU Governor restored to: $ORIGINAL_GOVERNOR"
    
    # Restore turbo boost (if it was available)
    if [ "$ORIGINAL_TURBO" != "N/A" ]; then
        echo $ORIGINAL_TURBO | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo > /dev/null
        echo "  Turbo setting restored to: $ORIGINAL_TURBO"
    fi
    
    # Restore ASLR
    echo $ORIGINAL_RANDOMIZE | sudo tee /proc/sys/kernel/randomize_va_space > /dev/null
    echo "  ASLR restored to: $ORIGINAL_RANDOMIZE"
    
    # Re-enable swap
    sudo swapon -a 2>/dev/null
    echo "  Swap re-enabled"
    
    echo -e "${GREEN}All settings restored!${NC}"
}

# Set up trap to restore settings on exit (Ctrl+C, etc.)
trap restore_settings EXIT INT TERM

# Apply benchmark settings
echo -e "${BLUE}Applying benchmark settings...${NC}"

# Set performance governor
echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor > /dev/null
echo "  CPU Governor: performance"

# Disable turbo boost (Intel)
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo > /dev/null 2>&1
echo "  Turbo boost: disabled"

# Disable ASLR
echo 0 | sudo tee /proc/sys/kernel/randomize_va_space > /dev/null
echo "  ASLR: disabled"

# Disable swap
sudo swapoff -a 2>/dev/null
echo "  Swap: disabled"

# Clear caches
sync
echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null
echo "  Caches: cleared"

echo -e "${GREEN}Benchmark environment ready!${NC}"
echo ""

# Parse command line arguments
T_VALUE=24  # Default value
RUNS=10    # Default value
BENCHMARK_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--power)
            T_VALUE="$2"
            shift 2
            ;;
        -r|--runs)
            RUNS="$2"
            shift 2
            ;;
        *)
            # Pass through any other arguments to the benchmark
            BENCHMARK_ARGS="$BENCHMARK_ARGS $1"
            shift
            ;;
    esac
done

# Run the actual benchmark
echo -e "${BLUE}Running benchmark with -t $T_VALUE -r $RUNS...${NC}"
sudo nice -n -20 taskset -c 0 ./bench1d_forward_fft bench -t $T_VALUE -r $RUNS $BENCHMARK_ARGS

echo ""
echo -e "${GREEN}Benchmark complete!${NC}"

# The trap will automatically call restore_settings() here