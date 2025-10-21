#ifndef CLI_PARSER_H
#define CLI_PARSER_H

#include <stdbool.h>

/**
 * @brief Defines the type of an argument's value.
 */
typedef enum {
  ARG_TYPE_BOOL,  // A flag that is either present (true) or not (false)
  ARG_TYPE_INT,   // An integer value (e.g., -n 10)
  ARG_TYPE_STRING // A string value (e.g., -f path/to/file)
} ArgType;

/**
 * @brief A structure to declaratively define a single command-line option.
 */
typedef struct {
  char short_name;         // The short name character (e.g., 'f' for -f)
  const char *long_name;   // The long name string (e.g., "file" for --file)
  const char *description; // The help text for this option
  ArgType type;            // The data type of the option's value
  void *value_ptr;  // Pointer to the variable where the parsed value will be
                    // stored
  bool is_required; // Whether this argument must be provided by the user
} CliOption;

/**
 * @brief Parses command-line arguments based on a declarative definition.
 *
 * This is the main function of the library. It processes argc/argv, populates
 * the variables pointed to by `value_ptr` in the options array, and handles
 * errors and help requests.
 *
 * @param argc The argument count from main().
 * @param argv The argument vector from main().
 * @param options An array of CliOption structs defining all possible arguments.
 * @param num_options The number of elements in the `options` array.
 * @param app_description A brief description of the application, printed at the
 * top of the help message.
 * @return 0 on successful parsing, -1 on error (e.g., missing required arg), 1
 * if help was printed.
 *
 * @note For ARG_TYPE_STRING, this function allocates memory for the string
 * using strdup(). The caller is responsible for freeing this memory.
 */
int cli_parse(int argc, char **argv, const CliOption options[], int num_options,
              const char *app_description);

/**
 * @brief Prints a formatted help message based on the option definitions.
 *
 * This function is automatically called by cli_parse() when -h or --help is
 * encountered, but can also be called manually.
 *
 * @param options An array of CliOption structs defining all possible arguments.
 * @param num_options The number of elements in the `options` array.
 * @param app_description A brief description of the application.
 */
void cli_print_help(const CliOption options[], int num_options,
                    const char *app_description);

#endif // CLI_PARSER_H