#include <iostream>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

#include <armadillo>

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkXMLPolyDataWriter.h"

//
// Types
//
typedef std::set < size_t >       VertexSetType;
typedef VertexSetType::iterator VertexSetItType;
typedef std::vector < VertexSetType > GraphType;
typedef arma::SpMat < double >    SparseMatType;
typedef arma::Mat   < double >     DenseMatType;
typedef arma::Col   < double >       VectorType;

//
// Main
//
int main( int argc, char ** argv ) {
  size_t num_vertices = atoi( argv[1] );
  size_t min_degree   = atoi( argv[2] );
  size_t offset       = atoi( argv[3] );
  char * out_file     =       argv[4];
  
  //
  // Make Random Graph
  //
  std::default_random_engine generator;
  std::uniform_int_distribution<size_t>
    distribution( 0, num_vertices - 1 );
  GraphType G = GraphType( num_vertices );
  for ( size_t i = 0; i < num_vertices; ++i ) {
    while ( G[i].size() < min_degree ) {
      size_t j = distribution(generator);
      if ( i != j ) {
        G[i].insert(j);
        G[j].insert(i);
      }
    }
  }

  //
  // Build Laplacian
  //
  SparseMatType L( num_vertices, num_vertices );
  for ( size_t i = 0; i < num_vertices; ++i ) {
    L( i, i ) = -G[i].size();
    for (
      VertexSetItType sit = G[i].begin();
      sit != G[i].end();
      ++sit
    ) {
      L( i, *sit ) = 1;
    }
  }

  //
  // Compute Eigen Values/Vectors
  //
  DenseMatType eigen_vectors;
  VectorType   eigen_values;
  eigs_sym( eigen_values, eigen_vectors, L, offset + 4, "sm" );

  //
  // Make output
  //
  vtkPoints * points = vtkPoints::New();
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints( num_vertices );
  for ( size_t i = 0; i < num_vertices; ++i ) {
    points->SetPoint( i,
      eigen_vectors( i, offset + 1 ),
      eigen_vectors( i, offset + 2 ),
      eigen_vectors( i, offset + 3 ));
  }
  vtkPolyData * graph_out = vtkPolyData::New();
  graph_out->Allocate();
  graph_out->SetPoints( points );
  vtkIdList * vtk_edge = vtkIdList::New();
  vtk_edge->SetNumberOfIds(2);
  for ( size_t i = 0; i < num_vertices; ++i ) {
    vtk_edge->SetId( 0, i );
    for (
      VertexSetItType sit = G[i].begin();
      sit != G[i].end();
      ++sit
    ) {
      size_t j = *sit;
      if ( i < j ) {
        vtk_edge->SetId( 1, j );
      }
      graph_out->InsertNextCell( VTK_LINE, vtk_edge );
    }
  }
  vtkXMLPolyDataWriter * writer = vtkXMLPolyDataWriter::New();
  writer->SetFileName( out_file );
  writer->SetInputData( graph_out );
  writer->Write();  
}
