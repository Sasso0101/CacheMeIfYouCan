// cli_parser.c

#define _GNU_SOURCE // For getopt_long
#include "cli_parser.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Internal helper to parse an integer and store it.
static int parse_int_arg(const char *optarg, int *dest) {
  char *endptr;
  long val = strtol(optarg, &endptr, 10);
  if (*endptr != '\0' || optarg == endptr) {
    fprintf(stderr, "Error: Invalid integer value '%s'\n", optarg);
    return -1;
  }
  *dest = (int)val;
  return 0;
}

void cli_print_help(const CliOption options[], int num_options,
                    const char *app_description) {
  printf("Usage: [options]\n");
  if (app_description) {
    printf("\n%s\n", app_description);
  }
  printf("\nOptions:\n");

  for (int i = 0; i < num_options; ++i) {
    const CliOption *opt = &options[i];
    char option_str[50];
    snprintf(option_str, sizeof(option_str), "-%c, --%s", opt->short_name,
             opt->long_name);

    const char *arg_type_str = "";
    if (opt->type == ARG_TYPE_INT)
      arg_type_str = " <int>";
    if (opt->type == ARG_TYPE_STRING)
      arg_type_str = " <string>";

    printf("  %-25s%s\t%s%s\n", option_str, arg_type_str, opt->description,
           opt->is_required ? " (required)" : "");
  }
  printf("  %-25s%s\t%s\n", "-h, --help", "",
         "Display this help message and exit");
}

int cli_parse(int argc, char **argv, const CliOption options[], int num_options,
              const char *app_description) {
  // 1. Construct the short options string for getopt_long
  // Example: "f:n:c" (f and n require arguments, c does not)
  char short_opts[num_options * 2 + 2]; // Max size needed
  int short_opts_idx = 0;
  for (int i = 0; i < num_options; ++i) {
    short_opts[short_opts_idx++] = options[i].short_name;
    if (options[i].type != ARG_TYPE_BOOL) {
      short_opts[short_opts_idx++] = ':';
    }
  }
  short_opts[short_opts_idx++] = 'h'; // Add 'h' for help
  short_opts[short_opts_idx] = '\0';

  // 2. Construct the long options array for getopt_long
  struct option long_opts[num_options + 2]; // +2 for help and null terminator
  for (int i = 0; i < num_options; ++i) {
    long_opts[i].name = options[i].long_name;
    long_opts[i].has_arg =
        (options[i].type != ARG_TYPE_BOOL) ? required_argument : no_argument;
    long_opts[i].flag = NULL;
    long_opts[i].val = options[i].short_name;
  }
  // Add help option
  long_opts[num_options] = (struct option){"help", no_argument, NULL, 'h'};
  // Null terminator for the array
  long_opts[num_options + 1] = (struct option){0, 0, 0, 0};

  // 3. Keep track of which options were seen
  bool *seen_options = calloc(num_options, sizeof(bool));
  if (!seen_options) {
    perror("calloc failed");
    return -1;
  }

  // 4. The parsing loop
  int opt_char;
  int long_index = 0;
  while ((opt_char = getopt_long(argc, argv, short_opts, long_opts,
                                 &long_index)) != -1) {
    if (opt_char == '?' || opt_char == ':') {
      free(seen_options);
      return -1; // getopt_long already printed an error
    }

    if (opt_char == 'h') {
      cli_print_help(options, num_options, app_description);
      free(seen_options);
      return 1; // Indicate that help was printed
    }

    // Find the matching option in our definitions array
    bool found = false;
    for (int i = 0; i < num_options; ++i) {
      if (options[i].short_name == opt_char) {
        seen_options[i] = true;
        found = true;
        const CliOption *opt_def = &options[i];
        switch (opt_def->type) {
        case ARG_TYPE_BOOL:
          *(bool *)opt_def->value_ptr = true;
          break;
        case ARG_TYPE_INT:
          if (parse_int_arg(optarg, (int *)opt_def->value_ptr) != 0) {
            free(seen_options);
            return -1;
          }
          break;
        case ARG_TYPE_STRING:
          // Caller is responsible for freeing this memory
          *(char **)opt_def->value_ptr = strdup(optarg);
          if (*(char **)opt_def->value_ptr == NULL) {
            perror("strdup failed");
            free(seen_options);
            return -1;
          }
          break;
        }
        break; // Found it, exit inner loop
      }
    }
    if (!found) {
      // Should not happen with well-formed options
      fprintf(stderr, "Warning: Parsed unknown option character '%c'.\n",
              opt_char);
    }
  }

  // 5. Check for missing required arguments
  int result = 0;
  for (int i = 0; i < num_options; ++i) {
    if (options[i].is_required && !seen_options[i]) {
      fprintf(stderr, "Error: Required argument '--%s' ('-%c') is missing.\n",
              options[i].long_name, options[i].short_name);
      result =
          -1; // Mark as error, but continue checking for other missing args
    }
  }

  free(seen_options);
  return result;
}