#include <iostream>

#include "CLI11.hpp"

#include "index_main.hpp"
#include "classify_main.hpp"
#include "version.hpp"

class MyFormatter : public CLI::Formatter {
public:
    std::string make_option_opts(const CLI::Option* opt) const override
    {
        std::stringstream out;

        if (opt->get_type_size() != 0) {
            if (!opt->get_type_name().empty())
                out << " " << get_label(opt->get_type_name());
            else if (opt->get_expected_min() > 1)
                out << " x " << opt->get_expected();

            if (opt->get_required())
                out << " " << get_label("[required]");
        }
        if (!opt->get_envname().empty())
            out << " (" << get_label("Env") << ":" << opt->get_envname() << ")";
        if (!opt->get_needs().empty()) {
            out << " " << get_label("Needs") << ":";
            for (const CLI::Option* op : opt->get_needs())
                out << " " << op->get_name();
        }
        if (!opt->get_excludes().empty()) {
            out << " " << get_label("Excludes") << ":";
            for (const CLI::Option* op : opt->get_excludes())
                out << " " << op->get_name();
        }
        return out.str();
    }

    std::string make_option_desc(const CLI::Option* opt) const override
    {
        std::stringstream out;
        out << opt->get_description();
        if (!opt->get_default_str().empty()) {
            out << " [default: " << opt->get_default_str() << "]";
        }
        return out.str();
    }
};

int main(int argc, char* argv[])
{
    CLI::App app { "Charon: Categorize reads into a small number of classes." };
    app.formatter(std::make_shared<MyFormatter>());
    app.add_flag_callback("-V,--version", []() {
        std::cout << "charon version " << CHARON_VERSION << std::endl;
        throw(CLI::Success {});
    });
    setup_index_subcommand(app);
    setup_classify_subcommand(app);


    app.require_subcommand();

    CLI11_PARSE(app, argc, argv);

    return 0;
}