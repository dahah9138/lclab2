#include <lclab2.h>
using AppSolver = LC::FrankOseen::Electric::RBFFDSolver;
using Dataset = AppSolver::Dataset;

class Sandbox : public LC::ConsoleApplication
{
public:
    explicit Sandbox(LC::Arguments arguments);
    void Run() override;
    ~Sandbox();

private:
};

Sandbox::Sandbox(LC::Arguments arguments) : LC::ConsoleApplication{arguments.AddRequired('f', "File").AddOptional('d', "Directory")} {
    _solver = std::make_unique<AppSolver>();
    Dataset* data = (Dataset*)(_solver->GetDataPtr());

    /*
    Q = 1 : Stable (dop = 3, V = 2)
    Q = 2 : Stable (dop = 5, V = 2.7-2.8)
    Q = 3 : Stable (dop = 7, V = 3.13)
    */
    int Q = 1;
    LC::scalar voltage = 2;
    LC::scalar dop = 2. * Q + 1.;

    int knn = 50;
    LC::scalar rNode = 0.14285714285;
    data->npp = 30.;
    LC::scalar qRatio = dop / (2. * Q + 1.);
    LC::scalar active_rad = Q + 0.5;

    (*data)
        .ElectricConstants(LC::FrankOseen::LC_TYPE::_5CB)
        .VoltageConfiguration(LC::Math::VoltageZ(0.0, voltage, dop))
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        .Rate(-.50)
        .Cell(dop, dop, dop)
        .Boundaries(0, 0, 0)
        .Neighbors(knn)
        .ExclusionRadius(LC::Math::ZLine(rNode, 0.25, 0.5, 1.5))
        .DirectorConfiguration(LC::Math::Heliknoton(Q, { dop, dop, dop }, 1. / qRatio));
		
    _solver->Init();

    LC_INFO("Created client application!");
}

void Sandbox::Run() {
    // Run the simulation
    bool GPU_ON = true;
    _solver->Relax(1000, GPU_ON);

    // Save the results
    std::string file = _arguments.Get('f');
    if (file.empty()) file = "default-toron.lmt";
    std::string directory = _arguments.Get('d');
    if (directory.empty()) directory = std::string(LCLAB2_ROOT_PATH) + "/data";
    _header.writeFile = directory + "/" + file;
    _solver->Export(_header);

    LC_INFO("Saved to file {0}", _header.writeFile.c_str());
}

Sandbox::~Sandbox() {
    LC_INFO("Destroying client application.");
}

LC::ConsoleApplication* LC::createConsoleApplication(int argc, char** argv) {
    return new Sandbox{ LC::Arguments{argc, argv} };
}