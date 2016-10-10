
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
//        qk::example::maxwell::gaussian_pulse("/home/seamill/local_storage/quikin/maxwell/gaussian_pulse/data");
//        qk::example::maxwell::wave_1D("/home/seamill/local_storage/quikin/maxwell/wave_1D/data");
//        qk::example::euler::sod_shock_tube("/home/seamill/local_storage/quikin/euler/sod_shock_tube/data");
//        qk::example::plasma::two_fluid_brio_wu("/home/seamill/local_storage/quikin/plasma/brio_wu/data");
//        qk::example::plasma::harris_current_sheet("/home/seamill/local_storage/quikin/plasma/harris_current_sheet/data");
        qk::example::plasma::gem_challenge("/home/seamill/local_storage/quikin/plasma/gem_challenge/data");
    } catch (const qk::exception & qke) {
        std::cout << "*** quikin exception caught ***\n" << qke.what();
        exit(EXIT_FAILURE);
    }
}
