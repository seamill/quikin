
// QK includes
#include "example/maxwell/gaussian_pulse.h"
#include "example/maxwell/wave_1D.h"
#include "example/euler/sod_shock_tube.h"
#include "example/plasma/two_fluid_brio_wu.h"
#include "example/plasma/harris_current_sheet.h"
#include "example/plasma/gem_challenge.h"
#include "lib/exception.h"

// STL includes
#include <iostream>

int main()
{
    try {
//        qk::example::maxwell::gaussian_pulse();
//        qk::example::maxwell::wave_1D();
//        qk::example::euler::sod_shock_tube();
//        qk::example::plasma::two_fluid_brio_wu();
//        qk::example::plasma::harris_current_sheet();
        qk::example::plasma::gem_challenge();
    } catch (qk::exception & qke) {
        std::cout << "*** quikin exception caught ***\n" << qke.what();
        exit(EXIT_FAILURE);
    }
}
