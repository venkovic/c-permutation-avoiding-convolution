#!/bin/bash
# diagnostic_bench1d.sh - Deep diagnosis of performance issues

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}=== SYSTEM PERFORMANCE DIAGNOSTIC ===${NC}"

# Function to check CPU frequency
check_cpu_freq() {
    echo -e "${YELLOW}CPU Frequency Status:${NC}"
    if command -v cpupower &> /dev/null; then
        cpupower frequency-info | grep -E "(current CPU frequency|governor)"
    else
        cat /proc/cpuinfo | grep "cpu MHz" | head -4
    fi
    echo ""
}

# Function to check thermal throttling
check_thermal() {
    echo -e "${YELLOW}Thermal Status:${NC}"
    if [ -f /sys/class/thermal/thermal_zone0/temp ]; then
        for zone in /sys/class/thermal/thermal_zone*/temp; do
            temp=$(cat $zone)
            temp_c=$((temp / 1000))
            zone_name=$(basename $(dirname $zone))
            echo "  $zone_name: ${temp_c}°C"
        done
    fi
    
    # Check for thermal throttling in dmesg
    if dmesg | tail -100 | grep -i "thermal\|throttl" | tail -3; then
        echo -e "${RED}  ⚠️  Recent thermal events detected${NC}"
    fi
    echo ""
}

# Function to check memory
check_memory() {
    echo -e "${YELLOW}Memory Status:${NC}"
    free -h
    echo ""
    
    echo -e "${YELLOW}Memory Frequency:${NC}"
    if command -v dmidecode &> /dev/null; then
        sudo dmidecode --type memory | grep -E "Speed|Configured" | head -4
    else
        echo "  dmidecode not available"
    fi
    echo ""
}

# Function to check system load
check_load() {
    echo -e "${YELLOW}System Load:${NC}"
    uptime
    echo ""
    
    echo -e "${YELLOW}Top CPU processes:${NC}"
    ps aux --sort=-%cpu | head -6
    echo ""
}

# Function to check power/performance settings
check_power() {
    echo -e "${YELLOW}Power Management:${NC}"
    
    # Check power profile (Ubuntu/modern systems)
    if command -v powerprofilesctl &> /dev/null; then
        echo "  Power profile: $(powerprofilesctl get)"
    fi
    
    # Check TLP (laptop power management)
    if command -v tlp-stat &> /dev/null; then
        echo "  TLP status: $(tlp-stat -s | grep "TLP status")"
    fi
    
    # Check cpupower
    if command -v cpupower &> /dev/null; then
        echo "  CPU governor: $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)"
    fi
    echo ""
}

# Function to run micro-benchmark
run_micro_benchmark() {
    echo -e "${YELLOW}Quick CPU micro-benchmark:${NC}"
    
    # Simple floating point benchmark
    start_time=$(date +%s.%N)
    awk 'BEGIN { for(i=0; i<1000000; i++) { x = i * 3.14159 * i } print "done" }' > /dev/null
    end_time=$(date +%s.%N)
    duration=$(echo "$end_time - $start_time" | bc -l)
    echo "  Floating point test: ${duration}s"
    
    # Memory bandwidth test
    start_time=$(date +%s.%N)
    dd if=/dev/zero of=/dev/null bs=1M count=1000 2>/dev/null
    end_time=$(date +%s.%N)
    duration=$(echo "$end_time - $start_time" | bc -l)
    echo "  Memory bandwidth test: ${duration}s"
    echo ""
}

# Main diagnostic
echo "Timestamp: $(date)"
echo "Hostname: $(hostname)"
echo "Kernel: $(uname -r)"
echo ""

check_cpu_freq
check_thermal
check_memory
check_load
check_power
run_micro_benchmark

# Parse command line arguments for benchmark
T_VALUE=20  # Smaller default for quick testing
RUNS=10

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
        --no-bench)
            echo -e "${BLUE}Skipping benchmark (--no-bench specified)${NC}"
            exit 0
            ;;
        *)
            shift
            ;;
    esac
done

echo -e "${BLUE}=== RUNNING ACTUAL BENCHMARK ===${NC}"
echo "Parameters: -t $T_VALUE -r $RUNS"
echo ""

# Check if benchmark exists
if [ ! -f "./bench1d" ]; then
    echo -e "${RED}Error: ./bench1d not found${NC}"
    echo "Please update the script with the correct benchmark executable name"
    exit 1
fi

# Run benchmark with monitoring
echo -e "${YELLOW}Monitoring during benchmark:${NC}"

# Start background monitoring
(
    while true; do
        temp=$(cat /sys/class/thermal/thermal_zone0/temp 2>/dev/null | head -1)
        freq=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq 2>/dev/null)
        if [ -n "$temp" ] && [ -n "$freq" ]; then
            temp_c=$((temp / 1000))
            freq_mhz=$((freq / 1000))
            echo "$(date +%H:%M:%S) - Temp: ${temp_c}°C, Freq: ${freq_mhz}MHz"
        fi
        sleep 2
    done
) &
monitor_pid=$!

# Run the benchmark
timeout 60 sudo nice -n -20 taskset -c 0 ./bench1d bench -t $T_VALUE -r $RUNS

# Stop monitoring
kill $monitor_pid 2>/dev/null

echo ""
echo -e "${BLUE}=== POST-BENCHMARK STATUS ===${NC}"
check_cpu_freq
check_thermal

echo -e "${GREEN}Diagnostic complete!${NC}"
echo ""
echo -e "${YELLOW}If performance is still poor, consider:${NC}"
echo "1. Reboot and test immediately"
echo "2. Check if system is thermal throttling"
echo "3. Verify no background processes are interfering"
echo "4. Test with different compiler optimizations"
echo "5. Check if hyperthreading/SMT is affecting results"