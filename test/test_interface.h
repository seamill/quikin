#ifndef _qk_tests_H
#define _qk_tests_H

// QK includes
#include "lib/exception.h"

// STL includes
#include <string>
#include <sstream>


namespace qk
{
namespace test
{

class test_monitor
{
public:
    test_monitor()
    {

    }

    ~test_monitor()
    {

    }

    void assert(const bool passed, const std::string & assertion_name)
    {
        if(!passed){
            throw qk::exception(assertion_name);
        }
        _log << "        Assertion '"<<assertion_name<<"' has passed.\n";
    }

    const std::string results() const {return _log.str();}
    void clear() {_log.str(std::string());}

protected:

    std::stringstream _log;

};


#define QK_TEST_CALL(test_name, test_function)                                      \
    test_stream << "    Testing '" << test_name << "'...\n";                        \
    try{                                                                            \
        const bool test_status = test_function(tm);                                 \
        test_passed = test_passed && test_status;                                   \
        test_stream << tm.results(); tm.clear();                                    \
        test_stream << "    Test '" << test_name << "' has "                        \
            << std::string((test_status) ? "succeeded." : "failed.") << std::endl;  \
    } catch(const std::exception & exc) {                                           \
        test_passed = false;                                                        \
        test_stream << tm.results(); tm.clear();                                    \
        test_stream << "    Test '"<<test_name<<"' Failed : " << exc.what() << "\n";\
    }

#define QK_TEST_CALL_GROUP_START(test_group_name)                       \
int                                                                     \
main(){                                                                 \
    std::ostream & test_stream = std::cout;                             \
    qk::test::test_monitor tm;                                          \
    bool test_passed = true;                                            \
    test_stream << "Running Test Group '"<<test_group_name<<"'...\n";

#define QK_TEST_CALL_GROUP_END(test_group_name)                         \
    test_stream << "Test Group '"<<test_group_name<<"' complete.\n";    \
    exit((test_passed) ? EXIT_SUCCESS : EXIT_FAILURE);                  \
}








}
}

#endif // _qk_variable_variable_H
