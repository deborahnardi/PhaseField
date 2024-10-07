#include "../headers/FEM.h"

void FEM::deleteFromString(std::string &fullStr, std::string removeStr)
{
    size_t found = fullStr.find(removeStr);
    if (found != std::string::npos)
        fullStr.erase(found, removeStr.length());
}

void FEM::writeInHDF5(hid_t &file, std::fstream &output_v, herr_t &status, hid_t &dataset, hid_t &dataspace, std::string AttributeName, std::string AttributeType, double *valueVector, hsize_t valueVectorDims[], std::string s1)
{
    output_v << "    <Attribute Name=\"" << AttributeName << "\" Center=\"Node\" AttributeType=\"" << AttributeType << "\" >" << std::endl;
    if (AttributeType == "Tensor")
        output_v << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"" << numNodes << " 9\">" << std::endl;
    else if (AttributeType == "Vector")
        output_v << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"" << numNodes << " 3\">" << std::endl;
    else if (AttributeType == "Scalar")
        output_v << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"" << numNodes << " 1\">" << std::endl;

    deleteFromString(AttributeName, " ");
    std::transform(AttributeName.begin(), AttributeName.end(), AttributeName.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    deleteFromString(s1, resultsPath + "results/");
    output_v << "        " << s1 << ":/" << AttributeName << std::endl;

    char name[AttributeName.size() + 2];
    std::snprintf(name, sizeof(name), "/%s", AttributeName.c_str());

    dataspace = H5Screate_simple(2, valueVectorDims, NULL);
    dataset = H5Dcreate2(file, name, H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &valueVector[0]);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    output_v << "      </DataItem>" << std::endl
             << "    </Attribute>" << std::endl;
}

void FEM::showResults()
{
    PetscPrintf(PETSC_COMM_WORLD, "Exporting data to Paraview...\n");
    int numElem = elements.size();
    int *connectivity = new int[numElNodes * numElem](); // () initializes the array with zeros
    double *tensor = new double[9 * numNodes]();
    double *vector = new double[3 * numNodes]();
    double *scalar = new double[numNodes]();
    /*
        .xdmf file is a common file format for visualization and data exchange in scientific computing, xdmf stands for eXtensible Data Model and Format;
        std::stream is a base class for all input/output streams, differently than std::ofstream, that is a base class for output streams only;
        std::ios_base::out is a flag that specifies output mode, it is used to open a file for writing;
    */

    std::string result;
    std::string r0 = resultsPath + "results/FEM_" + name + ".xdmf";
    std::fstream output_v(r0.c_str(), std::ios_base::out);

    std::string r1 = resultsPath + "results/hdf5/FEM_" + name + ".h5";
    std::fstream output_h5(r1.c_str(), std::ios_base::out);

    std::string r2 = resultsPath + "results/hdf5/FEM_" + name + "_connectivity.h5";

    std::string topology;
    if (elemDim == 1)
        topology = "Polyline";
    else if (elemDim == 2)
        topology = "Triangle";

    /*
        hid_t: it is a type defined in the HDF5 library, it is used to represent an abstract HDF5 object, it is an integer type;
        herr_t: it is a type defined in the HDF5 library, it is used to represent an error code, it is an integer type;
        hsize_t: it is a type defined in the HDF5 library, it is used to represent the size of an object, it is an integer type;
    */
    hid_t file, file2, dataset, dataspace;
    herr_t error;
    hsize_t tensorDims[2] = {numNodes, 9};               // Dimensions of the tensor dataset is 2D, numNodes x 9
    hsize_t vectorDims[2] = {numNodes, 3};               // Dimensions of the vector dataset is 2D, numNodes x 3
    hsize_t scalarDims[1] = {numNodes};                  // Dimensions of the scalar dataset is 1D, numNodes
    hsize_t connectivityDims[2] = {numElem, numElNodes}; // Dimensions of the connectivity dataset is 2D, numElem x numElNodes

    file = H5Fcreate(r1.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file2 = H5Fcreate(r2.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
        Writing the xdmf file
    */

    output_v << "<?xml version=\"1.0\"?>" << std::endl
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl
             << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\" >" << std::endl
             << "<Domain>" << std::endl
             << "  <Grid>" << std::endl
             << "    <Topology TopologyType=\"" << topology << "\" NumberOfElements=\"" << numElem << "\" >" << std::endl
             << "      <DataItem Format=    \"HDF\" NumberType=\"int\" Dimensions=\"" << numElem << " " << numElNodes << "\" >" << std::endl;

    /*
        Writing the connectivity dataset
    */

    output_v << "        hdf5/FEM_" + name + "_connectivity.h5:/connec" << std::endl;
    dataspace = H5Screate_simple(2, connectivityDims, NULL);
    dataset = H5Dcreate2(file2, "/connec", H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    error = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connectivity);
    error = H5Dclose(dataset);
    error = H5Sclose(dataspace);
    error = H5Fclose(file2);

    output_v << "      </DataItem>" << std::endl
             << "    </Topology>" << std::endl
             << "    <Geometry GeometryType=\"XYZ\">" << std::endl
             << "      <DataItem Format=\"HDF\" NumberType=\"double\" Dimensions=\"" << numNodes << " 3\">" << std::endl;

    /*
        Writing coordinates dataset
    */

    output_v << "        hdf5/FEM_" + name + "_coordinates.h5:/coord" << std::endl;

    for (auto n : discritizedNodes)
    {
        int index = n->getIndex();
        for (int i = 0; i < 3; i++)
            vector[3 * index + i] = n->getInitialCoordinates()[i];
    }
    dataspace = H5Screate_simple(2, vectorDims, NULL);
    dataset = H5Dcreate2(file, "/coords", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    error = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vector); // &vector[0] would also be correct
    error = H5Dclose(dataset);
    error = H5Sclose(dataspace);

    output_v << "      </DataItem>" << std::endl
             << "    </Geometry>" << std::endl;

    for (auto n : discritizedNodes)
    {
        int index = n->getIndex();
        for (int i = 0; i < 3; i++)
            vector[3 * index + i] = n->getInitialCoordinates()[i] - n->getFinalDisplacement()[i];
    }
    writeInHDF5(file, output_v, error, dataset, dataspace, "Displacement", "Vector", vector, vectorDims, r1);

    output_v << "  </Grid>" << std::endl
             << "</Domain>" << std::endl
             << "</Xdmf>" << std::endl;

    delete[] connectivity;
    delete[] tensor;
    delete[] vector;
    delete[] scalar;

    PetscPrintf(PETSC_COMM_WORLD, "Data exported successfully!\n");
}