#include <lclab2.h>
using AppSolver = LC::FrankOseen::ElasticOnly::FOFDSolver;
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

    (*data).Voxels(60, 60, 60)
        .Boundaries(1, 1, 0)
        .Cell(1, 1, 1)
        .ElasticConstants(LC::FrankOseen::ElasticConstants::_5CB())
        //.Configuration(Dataset::Heliknoton(1, translations, 0.2, 3));
        .Configuration(Dataset::Toron());
    //.Configuration(custom_cfg);
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