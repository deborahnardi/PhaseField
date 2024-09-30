void SolidDomain::exportToParaview(const int& step)
{
	std::stringstream text;
	text << "results/" << "solidOutput" << step << ".vtu";
	std::ofstream file(text.str());

	unsigned int numberOfElements = 0;
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}

	bool mixed = false;
	if (materials_[0]->getType() == MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID ||
		materials_[0]->getType() == MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID  ) mixed = true;

	//header
	file << "<?xml version=\"1.0\"?>" << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">" << "\n"
		 << "  <UnstructuredGrid>" << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << numberOfElements
         << "\">" << "\n";

	//nodal coordinates
	file << "    <Points>" << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">" << "\n";
    for (Node*& n : nodes_)
	{
		std::vector<DegreeOfFreedom*> dofs = n->getDegreesOfFreedom();
		file << dofs[0]->getCurrentValue() << " " << dofs[1]->getCurrentValue() << " " << ((dimension_ == 3)? dofs[2]->getCurrentValue() : 0.0) << "\n";
	}
    file << "      </DataArray>" << "\n"
         << "    </Points>" << "\n";

	//element connectivity
	file << "    <Cells>" << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">" << "\n";
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (unsigned int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n";

	//offsets
	file << "      <DataArray type=\"Int32\""
		 << " Name=\"offsets\" format=\"ascii\">" << "\n";
	int aux = 0;
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n";

	//elements type
	file << "      <DataArray type=\"UInt8\" Name=\"types\" "
		 << "format=\"ascii\">" << "\n";
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n"
		 << "    </Cells>" << "\n";
	
	//nodal results
	file << "    <PointData>" <<"\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ <<"\" "
		 << "Name=\"Displacement\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentValue() - dofs[i]->getInitialValue() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ << "\" "
		 << "Name=\"Velocity\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentFirstTimeDerivative() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ << "\" "
		 << "Name=\"Acceleration\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentSecondTimeDerivative() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_*(dimension_+1)/2 << "\" "
		 << "Name=\"CauchyStress\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		double* cauchyStress = n->getCauchyStress();
		for (int i = 0; i < dimension_*(dimension_+1)/2; i++)
			file << cauchyStress[i] << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	if (mixed)
	{
		file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 << "\" "
			<< "Name=\"Pressure\" format=\"ascii\">" << "\n";
		for (Node*& n: nodes_)
		{
			file << n->getDegreeOfFreedom(dimension_)->getCurrentValue() << "\n";
		}
		file << "      </DataArray> " << "\n";
	}
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 <<"\" "
		 << "Name=\"PermutedIndex\" format=\"ascii\">" << "\n";
	for (Node*& node : nodes_)
	{
		file << node->getPermutedIndex() << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "    </PointData>" << "\n";

	//elemental results
	file << "    <CellData>" << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 <<"\" "
		 << "Name=\"Rank\" format=\"ascii\">" << "\n";
	for (Element* const& el : elements_)
	{
		file << el->getRank() << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "    </CellData>" << "\n";

	//footnote
	file << "  </Piece>" << "\n"
		<< "  </UnstructuredGrid>" << "\n"
		<< "</VTKFile>" << "\n";
	file.close();

	// std::stringstream text2;
	// text2 << "error.dat";
	// std::ofstream file2;
	// if (timestep == 0) file2.open(text2.str());
	// else file2.open(text2.str(), std::ofstream::out | std::ofstream::app);
	// file2.precision(16);

	// file2 << timestep << "\n";
	// for (Node* n : nodes_)
	// {
	// 	file2 << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << "\n";
	// }

	// std::stringstream text2;
	// text2 << "error.dat";
	// std::ifstream file2;
	// file2.open(text2.str());
	// std::string line;
	// vector<double> reference_position(2 * nodes_.size());
	// for (int i = 0; i < 4000; i++)
	// {
	// 	std::getline(file2, line);
	// 	if (i == 16 * timestep + 15)
	// 	{
	// 		for (int j = 0; j < nodes_.size(); j++)
	// 		{
	// 			file2 >> reference_position(2 * j) >> reference_position(2 * j + 1);
	// 			std::getline(file2, line);
	// 		}
	// 		break;
	// 	}
	// 	else
	// 	{
	// 		for (int j = 0; j < nodes_.size(); j++)
	// 		{
	// 			std::getline(file2, line);
	// 		}
	// 	}
	// }

	// // Local error
	// std::stringstream text3;
	// text3 << "localerror-position" << deltat_ << ".dat";
	// std::ofstream file3;
	// file3.open(text3.str(), std::ofstream::out | std::ofstream::app);
	// double error_position = sqrt((nodes_[1]->getCurrentCoordinate()(0) - reference_position(2)) * (nodes_[1]->getCurrentCoordinate()(0) - reference_position(2)));

	// file3 << deltat_ * (timestep + 1.0) << " " << error_position << "\n";

	// // Norm L2
	// std::stringstream text4;
	// text4 << "globalerror-position" << deltat_ << ".dat";

	// std::ofstream file4;
	// file4.open(text4.str(), std::ofstream::out | std::ofstream::app);
	// error_position = 0.0;
	// for (Element* el : elements_)
	// {
	// 	vector<double> local_reference_position(2*el->getNodes().size());
	// 	for (int i = 0; i < el->getNodes().size(); i++)
	// 	{
	// 		int index = el->getNode(i)->getIndex();
	// 		local_reference_position(2*i) = reference_position(2*index);
	// 		local_reference_position(2*i+1) = reference_position(2*index+1);
	// 	}
	// 	error_position += el->domainIntegration(local_reference_position);
	// }
	// error_position = sqrt(error_position);
	// file4 << deltat_ * (timestep + 1.0) << " " << error_position << "\n";	
}