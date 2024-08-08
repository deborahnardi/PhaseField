Geometry *geo1 = new Geometry("rock_mesh1");


std::vector<Inclusion *> inclusions;
std::vector<MeshFactor *> factors;

inclusions.push_back(geo1->addInclusion(0.05, 0.2, 30., 0.25, 0.25));
inclusions.push_back(geo1->addInclusion(0.05, 0.2, 60., 0.75, 0.75));
inclusions.push_back(geo1->addInclusion(0.10, 0.10, 90., 0.5, 0.5));
inclusions.push_back(geo1->addInclusion(0.10, 0.10, 45., 0.25, 0.75));

factors.push_back(geo1->addMeshFactor(0.1, 1.0, 1.2, 1e-4, 1e-2));

geo1->setEdgeLength(1000.);
geo1->setAlgorithm(DELAUNAY);

std::cout << "Inclusions: " << inclusions.size() << std::endl;
std::cout << "Algorithm: " << geo1->getAlgorithm() << std::endl;
std::cout << "Edge length: " << geo1->getEdgeLength() << std::endl;
std::cout << "Mesh factors: " << factors.size() << std::endl;