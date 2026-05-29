/**
 * @file test_colors.h
 * @brief ANSI color escape codes for test output — shared by all test suites
 *
 * Usage:
 *   printf("%sFAILED%s: reason\n", COLOR_RED, COLOR_RESET);
 */

#pragma once

#define COLOR_RESET   "\033[0m"
#define COLOR_RED     "\033[1;31m"
#define COLOR_GREEN   "\033[1;32m"
#define COLOR_YELLOW  "\033[1;33m"
#define COLOR_CYAN    "\033[1;36m"
#define COLOR_BOLD    "\033[1m"
