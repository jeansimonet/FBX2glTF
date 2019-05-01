/**
 * Copyright (c) 2019-present, Adobe, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#include "Raw2Fbx.hpp"
#include <unordered_set>

static double scaleFactor = 100.0; // Technically we compute it in Raw2Fbx, but it will also come out to 100 :)

class Raw2FbxConverter
{
public:
  Raw2FbxConverter(const RawModel& rawModel, FbxScene* pFbxScene);
  bool WriteNodeHierarchy();

private:
  bool CreateNodeHierarchy();
  bool CreateTextures();
  bool CreateMaterials();
  FbxMesh* CreateBaseMesh(int surfaceIndex, std::vector<int>& surfaceVertIndices);
  bool CreateSkeleton(const RawSurface& surface);
  bool CreateSkin(const RawSurface& surface, const std::vector<int>& surfaceVertIndices);
  bool CreateMorphTargets(const RawSurface& surface, const std::vector<int>& surfaceVertIndices);
  bool CreateMeshes();
  bool AssignMeshesAndMaterialsToNodes();
  bool CreateAnimations();
  bool CreateLights();
  bool CreateCameras();

private:
  const RawModel& raw;
  FbxScene* pScene;

  // We need a small map of surface id to node(s) that use that surface
  // This is initialized in the constructor
  std::multimap<long, long> surfaceIdToNodeIdsMap;
  std::map<long, FbxNode*> idToFbxNodeMap;
  std::vector<FbxFileTexture*> fbxTextures;
  std::vector<FbxSurfaceMaterial*> fbxSurfaceMaterials;
  std::map<long, FbxMesh*> surfaceIdToFbxMeshMap;
  std::map<long, int> surfaceIdToMatIndexMap;
  int attributes;
};


Raw2FbxConverter::Raw2FbxConverter(const RawModel& rawModel, FbxScene* pFbxScene)
: raw(rawModel)
, pScene(pFbxScene)
, attributes(rawModel.GetVertexAttributes())
{
  // Initialize our surface lookup map
  for (int i = 0; i < raw.GetNodeCount(); ++i) {
    // For any node that has a surface Id, add a mesh to the matching fbx node
    auto& node = raw.GetNode(i);

    // Does the node have a surface
    if (node.surfaceId >= 0) {
        // Add/Append a new entry
        surfaceIdToNodeIdsMap.insert( {node.surfaceId, node.id } );
    }
  }
}

bool Raw2FbxConverter::CreateNodeHierarchy() {
  // Create all the matching fbx node of the raw model
  for (int i = 0; i < raw.GetNodeCount(); ++i) {

      auto& node = raw.GetNode(i);

      // Create the new node
      FbxNode* fbxNode = FbxNode::Create(pScene, node.name.c_str());

      // Add it to the map for later parent lookup
      idToFbxNodeMap[node.id] = fbxNode;

      // Set transform inheritance to the only kind we support right now
      fbxNode->SetTransformationInheritType(FbxTransform::eInheritRSrs);

      // Set the node local transform
      FbxAMatrix mat;
      mat.SetTQS(toFbxDouble3(node.translation), toFbxQuaternion(node.rotation), toFbxVector4Position(node.scale));
      fbxNode->LclTranslation.Set(mat.GetT());
      fbxNode->LclRotation.Set(mat.GetR());
      fbxNode->LclScaling.Set(mat.GetS());
  }

  // Set all the parent links
  for (int i = 0; i < raw.GetNodeCount(); ++i) {

      auto& node = raw.GetNode(i);

      // Find the parent node in the fbx, or root if none
      auto parentFbxNode = pScene->GetRootNode();
      auto parentFbxNodeIt = idToFbxNodeMap.find(node.parentId);
      if (parentFbxNodeIt != idToFbxNodeMap.end()) {
        parentFbxNode = parentFbxNodeIt->second;
      }

      // Find the node itself
      FbxNode* fbxNode = idToFbxNodeMap[node.id];
      if (!fbxNode) {
             fmt::printf("Error:: Cannot find node %s in created fbx tree.\n", node.name.c_str());
             return false;
       }

      // Add the child to its parent
      parentFbxNode->AddChild(fbxNode);
  }

  return true;
}

bool Raw2FbxConverter::CreateTextures() {
  // Create Textures
  for (int i = 0; i < raw.GetTextureCount(); ++i) {
    auto& rawTexture = raw.GetTexture(i);
    FbxFileTexture* pTexture = FbxFileTexture::Create(pScene, rawTexture.name.c_str());
    pTexture->SetFileName(rawTexture.fileLocation.c_str());

    // Set the type of texture
    switch (rawTexture.usage) {
      case RAW_TEXTURE_USAGE_NORMAL:
        pTexture->SetTextureUse(FbxTexture::eBumpNormalMap);
        break;
      case RAW_TEXTURE_USAGE_REFLECTION:
        pTexture->SetTextureUse(FbxTexture::eSphericalReflectionMap);
        break;
      default:
        pTexture->SetTextureUse(FbxTexture::eStandard);
        break;
    }

    pTexture->SetMappingType(FbxTexture::eUV);
    pTexture->SetMaterialUse(FbxFileTexture::eModelMaterial);

    fbxTextures.push_back(pTexture);
  }

  return true;
}

bool Raw2FbxConverter::CreateMaterials() {
  // Materials
  for (int i = 0; i < raw.GetMaterialCount(); ++i) {
    auto& rawMaterial = raw.GetMaterial(i);
    
    // Create the proper surface type
    FbxSurfaceMaterial* pFbxSurf = nullptr;
    switch (rawMaterial.info->shadingModel) {
      case RAW_SHADING_MODEL_CONSTANT:
        // TODO
        break;
      case RAW_SHADING_MODEL_LAMBERT:
        {
          // Create the surface
          FbxSurfaceLambert* pFbxSurfLambert = FbxSurfaceLambert::Create(pScene, rawMaterial.name.c_str());
          auto props = std::static_pointer_cast<RawTraditionalMatProps>(rawMaterial.info);
          
          // Set material textures and constants
          if (rawMaterial.textures[RAW_TEXTURE_USAGE_EMISSIVE] != -1) {
            pFbxSurfLambert->Emissive.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_EMISSIVE]]);
          } else {
            pFbxSurfLambert->Emissive = toFbxDouble3(props->emissiveFactor);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_AMBIENT] != -1) {
            pFbxSurfLambert->Ambient.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_AMBIENT]]);
          } else {
            pFbxSurfLambert->Ambient = toFbxDouble3(props->ambientFactor);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_DIFFUSE] != -1) {
            pFbxSurfLambert->Diffuse.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_DIFFUSE]]);
          } else {
            pFbxSurfLambert->Diffuse = toFbxVector4(props->diffuseFactor);
          }

          pFbxSurf = pFbxSurfLambert;
        }
        break;
      case RAW_SHADING_MODEL_BLINN:
        // TODO
        break;
      case RAW_SHADING_MODEL_PHONG:
        {
          // Create the surface
          FbxSurfacePhong* pFbxSurfPhong = FbxSurfacePhong::Create(pScene, rawMaterial.name.c_str());
          auto props = std::static_pointer_cast<RawTraditionalMatProps>(rawMaterial.info);
          
          // Set material textures and constants
          // Set material textures and constants
          if (rawMaterial.textures[RAW_TEXTURE_USAGE_EMISSIVE] != -1) {
            pFbxSurfPhong->Emissive.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_EMISSIVE]]);
          } else {
            pFbxSurfPhong->Emissive = toFbxDouble3(props->emissiveFactor);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_AMBIENT] != -1) {
            pFbxSurfPhong->Ambient.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_AMBIENT]]);
          } else {
            pFbxSurfPhong->Ambient = toFbxDouble3(props->ambientFactor);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_DIFFUSE] != -1) {
            pFbxSurfPhong->Diffuse.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_DIFFUSE]]);
          } else {
            pFbxSurfPhong->Diffuse = toFbxVector4(props->diffuseFactor);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_NORMAL] != -1) {
            pFbxSurfPhong->NormalMap.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_NORMAL]]);
          }

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_SPECULAR] != -1) {
            pFbxSurfPhong->Specular.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_SPECULAR]]);
          } else {
            pFbxSurfPhong->Specular = toFbxDouble3(props->specularFactor);;
          }
          pFbxSurfPhong->Shininess = props->shininess;

          if (rawMaterial.textures[RAW_TEXTURE_USAGE_REFLECTION] != -1) {
            pFbxSurfPhong->Reflection.ConnectSrcObject(fbxTextures[rawMaterial.textures[RAW_TEXTURE_USAGE_REFLECTION]]);
          }

          pFbxSurf = pFbxSurfPhong;
        }
        break;
      case RAW_SHADING_MODEL_PBR_MET_ROUGH:
        // TODO
        break;
      default:
        break;
    }

    // Add the material to the map so we can connect it to the meshes that use it!
    fbxSurfaceMaterials.push_back(pFbxSurf);
  }

  return true;
}

FbxMesh* Raw2FbxConverter::CreateBaseMesh(int surfaceIndex, std::vector<int>& surfaceVertIndices) {

  auto& surface = raw.GetSurface(surfaceIndex);

  FbxMesh* fbxMesh = FbxMesh::Create(pScene, surface.name.c_str());
  
  // We need to collect all the triangles and verts of this surface

  // We'll store the indices of the triangles and verts used by this surface
  std::vector<int> surfaceTriIndices;

  // And a mapping of the verts index in the raw model to the vert index in the collected array above
  // For this we iterate over all the triangles in the model and select the ones that match this surface.
  // With a lot of surfaces this will be very inefficient, and should instead split the process up
  // with more intermediate data and lookup.
  std::map<int, int> surfaceVertIndexToMeshVertIndex;
  for (int t = 0; t < raw.GetTriangleCount(); ++t)
  {
    auto& tri = raw.GetTriangle(t);
    if (tri.surfaceIndex == surfaceIndex) {
      surfaceTriIndices.push_back(t);

      // Add / Update the material index for this surface
      surfaceIdToMatIndexMap[surface.id] = tri.materialIndex;

      for (int v = 0; v < 3; ++v) {
        // Have we already seen this vert?
        auto&& prevVertIt = surfaceVertIndexToMeshVertIndex.find(tri.verts[v]);
        if (prevVertIt == surfaceVertIndexToMeshVertIndex.end()) {
          // We have not! add it! First Add the index lookup
          surfaceVertIndexToMeshVertIndex[tri.verts[v]] = surfaceVertIndices.size();

          // And add to the vert list
          surfaceVertIndices.push_back(tri.verts[v]);
        }
        // Else we've seen this vert already
      }
    }
    // Else this triangle doesn't interest us
  }

  // Let's start by adding all the verts to the mesh
  fbxMesh->InitControlPoints(surfaceVertIndices.size());

  // Add geometry elements based on the vertex attributes
  FbxVector4* controlPoints = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_POSITION) {
    controlPoints = fbxMesh->GetControlPoints();
  }

  FbxGeometryElementNormal* lGeometryElementNormal = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_NORMAL) {
    lGeometryElementNormal = fbxMesh->CreateElementNormal();
    lGeometryElementNormal->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementNormal->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  FbxGeometryElementTangent* lGeometryElementTangent = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_TANGENT) {
    lGeometryElementTangent = fbxMesh->CreateElementTangent();
    lGeometryElementTangent->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementTangent->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  FbxGeometryElementBinormal* lGeometryElementBinormal = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_BINORMAL) {
    lGeometryElementBinormal = fbxMesh->CreateElementBinormal();
    lGeometryElementBinormal->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementBinormal->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  FbxGeometryElementVertexColor* lGeometryElementVertexColor = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_COLOR) {
    lGeometryElementVertexColor = fbxMesh->CreateElementVertexColor();
    lGeometryElementVertexColor->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementVertexColor->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  FbxGeometryElementUV* lGeometryElementUV0 = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_UV0) {
    // Not sure that I am using the correct layer element type here...
    lGeometryElementUV0 = fbxMesh->CreateElementUV((surface.name + "UV0").c_str(), FbxLayerElement::eTextureDiffuse);
    lGeometryElementUV0->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementUV0->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  FbxGeometryElementUV* lGeometryElementUV1 = nullptr;
  if (attributes & RAW_VERTEX_ATTRIBUTE_UV1) {
    // Not sure that I am using the correct layer element type here...
    lGeometryElementUV1 = fbxMesh->CreateElementUV((surface.name + "UV1").c_str(), FbxLayerElement::eUV);
    lGeometryElementUV1->SetMappingMode(FbxGeometryElement::eByControlPoint);
    lGeometryElementUV1->SetReferenceMode(FbxGeometryElement::eDirect);
  }

  for (int v = 0; v < surfaceVertIndices.size(); ++v) {
    auto& vert = raw.GetVertex(surfaceVertIndices[v]);

    if (attributes & RAW_VERTEX_ATTRIBUTE_POSITION) {
      controlPoints[v] = toFbxVector4Position(vert.position);
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_NORMAL) {
      lGeometryElementNormal->GetDirectArray().Add(toFbxVector4Vector(vert.normal));
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_TANGENT) {
      lGeometryElementTangent->GetDirectArray().Add(toFbxVector4(vert.tangent));
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_BINORMAL) {
      lGeometryElementBinormal->GetDirectArray().Add(toFbxVector4Vector(vert.binormal));
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_COLOR) {
      lGeometryElementVertexColor->GetDirectArray().Add(toFbxColor(vert.color));
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_UV0) {
      lGeometryElementUV0->GetDirectArray().Add(toFbxVector2(vert.uv0));
    }

    if (attributes & RAW_VERTEX_ATTRIBUTE_UV1) {
      lGeometryElementUV1->GetDirectArray().Add(toFbxVector2(vert.uv1));
    }
  }

  // Create Polygons for each triangle of the surface
  for (int t = 0; t < surfaceTriIndices.size(); ++t) {
    auto& tri = raw.GetTriangle(surfaceTriIndices[t]);
    fbxMesh->BeginPolygon(-1, -1, -1, false);
    fbxMesh->AddPolygon(surfaceVertIndexToMeshVertIndex[tri.verts[0]]);
    fbxMesh->AddPolygon(surfaceVertIndexToMeshVertIndex[tri.verts[1]]);
    fbxMesh->AddPolygon(surfaceVertIndexToMeshVertIndex[tri.verts[2]]);
    fbxMesh->EndPolygon();
  }

  return fbxMesh;
}

bool Raw2FbxConverter::CreateSkeleton(const RawSurface& surface) {
  // Create the skeletons
  auto& rawSkeletonRootNode = raw.GetNode(raw.GetNodeById(surface.skeletonRootId));
  if (!rawSkeletonRootNode.isJoint) {
    return false;
  }

  FbxSkeleton* pSkeletonRootAttribute = FbxSkeleton::Create(pScene, rawSkeletonRootNode.name.c_str());
  pSkeletonRootAttribute->SetSkeletonType(FbxSkeleton::eLimb);
  FbxNode* pSkeletonRoot = idToFbxNodeMap[surface.skeletonRootId];
  pSkeletonRoot->SetNodeAttribute(pSkeletonRootAttribute);

  // Now iterate through the children and children of children

  // Define a lambda to do the work of adding a skeleton attribute to a bone node and
  // recursing through the hierarchy of bones
  std::function<void(int)> attachSkeletonAttribute = [&] (int skeletonNodeId) {
    auto& rawSkeletonNode = raw.GetNode(raw.GetNodeById(skeletonNodeId));
    if (!rawSkeletonNode.isJoint) {
      // Skip non-joints in the joint chain.
      return;
    }

    // Add the skeleton attribute to the node
    FbxSkeleton* pSkeletonNodeAttribute = FbxSkeleton::Create(pScene, rawSkeletonNode.name.c_str());
    // SDK seems to indicate that the type should be eEffector for the last node in a tree, but
    // the examples don't follow that, so not sure if it is really necessary
    pSkeletonNodeAttribute->SetSkeletonType(FbxSkeleton::eLimb); 
    FbxNode* pSkeletonNode = idToFbxNodeMap[skeletonNodeId];
    pSkeletonNode->SetNodeAttribute(pSkeletonNodeAttribute);

    // Recurse!
    for (int c = 0; c < rawSkeletonNode.childIds.size(); ++c) {
      attachSkeletonAttribute(rawSkeletonNode.childIds[c]);
    }
  };

  // Start by doing the work on children of the root
  for (auto childId : rawSkeletonRootNode.childIds) {
    attachSkeletonAttribute(childId);
  }
  return true;
}

bool Raw2FbxConverter::CreateSkin(const RawSurface& surface, const std::vector<int>& surfaceVertIndices) {

  // Find the (first/only?) node that corresponds to the surface
  auto rawPatchNodeIdIt = surfaceIdToNodeIdsMap.find(surface.id);
  assert(rawPatchNodeIdIt != surfaceIdToNodeIdsMap.end());
  long rawPatchNodeId = rawPatchNodeIdIt->second;
  auto& rawPatchNode = raw.GetNode(raw.GetNodeById(rawPatchNodeId));

  FbxAMatrix transformMatrix;
  transformMatrix.SetTQS(
    toFbxVector4Vector(rawPatchNode.translation),
    toFbxQuaternion(rawPatchNode.rotation),
    toFbxVector4Vector(rawPatchNode.scale));

  // Build a map of verts affected by each bone
  std::map<long, FbxCluster*> rawBoneNodeIdToClusterMap;

  // To do this, we iterate all the verts of the surface,
  // adding the vert to the proper cluster as we go.
  for (int v = 0; v < surfaceVertIndices.size(); ++v) {
    auto& vert = raw.GetVertex(surfaceVertIndices[v]);

    // Up to 4 joints per vert
    for (int j = 0; j < 4; ++j) {
      int jointIndex = vert.jointIndices[j];
      if (jointIndex >= 0) {
        float jointWeight = vert.jointWeights[j];
        if (jointWeight > 0.0f) {
          // From the index, we want the node id, so we can look it up
          long skeletonNodeId = surface.jointIds[jointIndex];

          // Find the cluster for this joint
          FbxCluster* pCluster = nullptr;
          auto clusterIt = rawBoneNodeIdToClusterMap.find(skeletonNodeId);
          if (clusterIt == rawBoneNodeIdToClusterMap.end()) {
            // There isn't one already, create a cluster for this bone
            auto& rawSkeletonNode = raw.GetNode(raw.GetNodeById(skeletonNodeId));
            FbxNode* pSkeletonNode = idToFbxNodeMap[skeletonNodeId];

            pCluster = FbxCluster::Create(pScene, rawSkeletonNode.name.c_str());
            pCluster->SetLink(pSkeletonNode);
            pCluster->SetLinkMode(FbxCluster::eNormalize);

            // Set the transforms
            pCluster->SetTransformMatrix(transformMatrix);
            pCluster->SetTransformLinkMatrix(transformMatrix * toFbxAMatrix(surface.inverseBindMatrices[jointIndex].Inverse()));

            // Add the cluster to the map
            rawBoneNodeIdToClusterMap[skeletonNodeId] = pCluster;
          } else {
            // There was already a cluster for this bone
            pCluster = clusterIt->second;
          }

          assert(pCluster); // FIXME, Add error message
          pCluster->AddControlPointIndex(v, jointWeight);
        }
      }
    }
  }

  // Add the skin to the mesh, and the clusters to the skin!
  FbxSkin* pSkin = FbxSkin::Create(pScene, surface.name.c_str());
  pSkin->SetSkinningType(FbxSkin::eBlend);

  for (auto idAndCluster : rawBoneNodeIdToClusterMap) {
    pSkin->AddCluster(idAndCluster.second);
  }

  // Attach the skin deformer to the mesh
  FbxMesh* fbxMesh = surfaceIdToFbxMeshMap[surface.id];
  fbxMesh->AddDeformer(pSkin);

  return true;
}

bool Raw2FbxConverter::CreateMorphTargets(const RawSurface& surface, const std::vector<int>& surfaceVertIndices) {
  // Blend channels. This is kind of like skinning data, in that FBX works
  // the opposite way that RAW does. Raw stores the blends in each vertex,
  // while the FBX stores a modifier with a duplicate shape (i.e mesh) for
  // each channel.
  for (int c = 0; c < surface.blendChannels.size(); ++c) {
    // For each channel, we need to create a new shape and control points
    auto& channel = surface.blendChannels[c];

    // Creathe shape etc...
    FbxShape* pShape = FbxShape::Create(pScene, channel.name.c_str());
    pShape->InitControlPoints(surfaceVertIndices.size());

    // Add geometry elements based on the vertex attributes
    FbxVector4* controlPoints = nullptr;
    if (attributes & RAW_VERTEX_ATTRIBUTE_POSITION) {
      controlPoints = pShape->GetControlPoints();
    }

    FbxGeometryElementNormal* lGeometryElementNormal = nullptr;
    if (attributes & RAW_VERTEX_ATTRIBUTE_NORMAL) {
      lGeometryElementNormal = pShape->CreateElementNormal();
      lGeometryElementNormal->SetMappingMode(FbxGeometryElement::eByControlPoint);
      lGeometryElementNormal->SetReferenceMode(FbxGeometryElement::eDirect);
    }

    FbxGeometryElementTangent* lGeometryElementTangent = nullptr;
    if (attributes & RAW_VERTEX_ATTRIBUTE_TANGENT) {
      lGeometryElementTangent = pShape->CreateElementTangent();
      lGeometryElementTangent->SetMappingMode(FbxGeometryElement::eByControlPoint);
      lGeometryElementTangent->SetReferenceMode(FbxGeometryElement::eDirect);
    }

    for (int v = 0; v < surfaceVertIndices.size(); ++v) {
      auto& vert = raw.GetVertex(surfaceVertIndices[v]);

      if (attributes & RAW_VERTEX_ATTRIBUTE_POSITION) {
        controlPoints[v] = toFbxVector4Position(vert.blends[c].position + vert.position);
      }

      if (attributes & RAW_VERTEX_ATTRIBUTE_NORMAL) {
        lGeometryElementNormal->GetDirectArray().Add(toFbxVector4Vector(vert.blends[c].normal + vert.normal));
      }

      if (attributes & RAW_VERTEX_ATTRIBUTE_TANGENT) {
        lGeometryElementTangent->GetDirectArray().Add(toFbxVector4(vert.blends[c].tangent + vert.tangent));
      }
    }

    // Now associate this shape with a channel
    FbxBlendShape* pBlendShape = FbxBlendShape::Create(pScene, channel.name.c_str());
    FbxBlendShapeChannel* pBlendShapeChannel = FbxBlendShapeChannel::Create(pScene, channel.name.c_str());
    pBlendShape->AddBlendShapeChannel(pBlendShapeChannel);
    pBlendShapeChannel->AddTargetShape(pShape);

    // And finally the modifier to the mesh
    FbxMesh* fbxMesh = surfaceIdToFbxMeshMap[surface.id];
    fbxMesh->AddDeformer(pBlendShape);
  }
  return true;
}

bool Raw2FbxConverter::CreateMeshes() {
  for (int i = 0; i < raw.GetSurfaceCount(); ++i) {
    auto& surface = raw.GetSurface(i);

    // Create the base mesh for this surface
    std::vector<int> surfaceVertIndices;
    FbxMesh* fbxMesh = CreateBaseMesh(i, surfaceVertIndices);

    // Add the mesh to the map, so we can map it later to actual fbx nodes
    surfaceIdToFbxMeshMap[surface.id] = fbxMesh;

    // Skinning
    if (attributes & RAW_VERTEX_ATTRIBUTE_JOINT_INDICES || attributes & RAW_VERTEX_ATTRIBUTE_JOINT_WEIGHTS) {

      // First create the skeleton for this surface, if any
      CreateSkeleton(surface);

      CreateSkin(surface, surfaceVertIndices);
    }

    CreateMorphTargets(surface, surfaceVertIndices);
  }

  return true;
}

bool Raw2FbxConverter::AssignMeshesAndMaterialsToNodes() {
  // Now for each node that has a valid surface Id, associate the fbxMesh with the
  // correct fbxNode. This means a mesh can be associated with multiple nodes if desired.
  // I don't know if that's the right thing to do, but it seems valid based on the sdk.
  for (int i = 0; i < raw.GetNodeCount(); ++i) {
    // For any node that has a surface Id, add a mesh to the matching fbx node
    auto& node = raw.GetNode(i);
    // Does this node have a mesh? If so, grab it!
    FbxMesh* fbxMesh = surfaceIdToFbxMeshMap[node.surfaceId];
    if (fbxMesh) {
      // Grab matching fbxNode
      FbxNode* fbxNode = idToFbxNodeMap[node.id];
      if (!fbxNode) {
            fmt::printf("Error:: Cannot find node %s.\n", node.name.c_str());
            return false;
      }

      // Connect!!
      fbxNode->SetNodeAttribute(fbxMesh);

      // Assign materials
      // Right now, we only support a single material
      int firstMaterialIndex = surfaceIdToMatIndexMap[node.surfaceId];
      FbxSurfaceMaterial* fbxMat = fbxSurfaceMaterials[firstMaterialIndex];
      if (fbxMat) {
        FbxGeometryElementMaterial* pGeometryElementMaterial = fbxMesh->GetElementMaterial(0);
        if (!pGeometryElementMaterial) {
            pGeometryElementMaterial = fbxMesh->CreateElementMaterial();
        }

        // The material is mapped to the whole mesh
        pGeometryElementMaterial->SetMappingMode(FbxGeometryElement::eAllSame);

        // And the material is avalible in the Direct array
        pGeometryElementMaterial->SetReferenceMode(FbxGeometryElement::eDirect);
        fbxNode->AddMaterial(fbxMat);
      }
    }
  }
  return true;
}

bool Raw2FbxConverter::CreateAnimations() {
  // Animations
  for (int i = 0; i < raw.GetAnimationCount(); ++i) {
    // Get raw animation
    auto& rawAnim = raw.GetAnimation(i);
    
    // Create an anim stack in the fbx
    FbxAnimStack* pAnimStack = FbxAnimStack::Create(pScene, rawAnim.name.c_str());

    // Set the time range
    pAnimStack->LocalStart = rawAnim.times.front();
    pAnimStack->LocalStop = rawAnim.times.back();

    // Create one (and only) layer
    FbxAnimLayer* pAnimLayer = FbxAnimLayer::Create(pScene, "Base Layer");
    pAnimStack->AddMember(pAnimLayer);

    // Add a curve for each node and each property of the initial anim
    for (int c = 0; c < rawAnim.channels.size(); ++c) {
      auto& rawChannel = rawAnim.channels[c];

      // Grab the fbx bone this channel points to
      auto& rawNode = raw.GetNode(rawChannel.nodeIndex);
      FbxNode* pFbxNode = idToFbxNodeMap[rawNode.id];

      // TODO: Need to animate blend shapes weight...

      // Now add data to those!
      FbxTime keyTime;
      FbxAnimCurveKey key;

      // First the translations
      if (!rawChannel.translations.empty()) {
        // Tell the SDK to create the compound XYZ curve, we'll retrieve individual x y and z later.
        pFbxNode->LclTranslation.GetCurveNode(pAnimLayer, true);
        auto addTranslationCurve = [&] (int component, const char* channelId) {
          FbxAnimCurve* pCurve = pFbxNode->LclTranslation.GetCurve(pAnimLayer, channelId, true);
          pCurve->KeyModifyBegin();
          for (int t = 0; t < rawAnim.times.size(); ++t) {
            keyTime.SetSecondDouble(rawAnim.times[t]);
            key.Set(keyTime, toFbxDouble3(rawChannel.translations[t])[component]);
            pCurve->KeyAdd(keyTime, key);
          }
          pCurve->KeyModifyEnd();
        };
        addTranslationCurve(0, FBXSDK_CURVENODE_COMPONENT_X);
        addTranslationCurve(1, FBXSDK_CURVENODE_COMPONENT_Y);
        addTranslationCurve(2, FBXSDK_CURVENODE_COMPONENT_Z);
      }

      // Then rotations
      if (!rawChannel.rotations.empty()) {
        // Tell the SDK to create the compound XYZ curve, we'll retrieve individual x y and z later.
        pFbxNode->LclRotation.GetCurveNode(pAnimLayer, true);
        auto addRotationCurve = [&] (int component, const char* channelId) {
          FbxAnimCurve* pCurve = pFbxNode->LclRotation.GetCurve(pAnimLayer, channelId, true);
          pCurve->KeyModifyBegin();
          for (int t = 0; t < rawAnim.times.size(); ++t) {
            keyTime.SetSecondDouble(rawAnim.times[t]);

            // So ugly...
            FbxAMatrix helpMeFixThisRoitationPlease;
            helpMeFixThisRoitationPlease.SetTQS(FbxDouble3(0,0,0), toFbxQuaternion(rawChannel.rotations[t]), FbxVector4(1,1,1,1));
            key.Set(keyTime, helpMeFixThisRoitationPlease.GetR()[component]);
            pCurve->KeyAdd(keyTime, key);
          }
          pCurve->KeyModifyEnd();
        };
        addRotationCurve(0, FBXSDK_CURVENODE_COMPONENT_X);
        addRotationCurve(1, FBXSDK_CURVENODE_COMPONENT_Y);
        addRotationCurve(2, FBXSDK_CURVENODE_COMPONENT_Z);
      }

      // Finally Scaling
      if (!rawChannel.scales.empty()) {
        // Tell the SDK to create the compound XYZ curve, we'll retrieve individual x y and z later.
        pFbxNode->LclScaling.GetCurveNode(pAnimLayer, true);
        auto addScalingCurve = [&] (int component, const char* channelId) {
          FbxAnimCurve* pCurve = pFbxNode->LclScaling.GetCurve(pAnimLayer, channelId, true);
          pCurve->KeyModifyBegin();
          for (int t = 0; t < rawAnim.times.size(); ++t) {
            keyTime.SetSecondDouble(rawAnim.times[t]);
            key.Set(keyTime, toFbxVector4Position(rawChannel.scales[t])[component]);
            pCurve->KeyAdd(keyTime, key);
          }
          pCurve->KeyModifyEnd();
        };
        addScalingCurve(0, FBXSDK_CURVENODE_COMPONENT_X);
        addScalingCurve(1, FBXSDK_CURVENODE_COMPONENT_Y);
        addScalingCurve(2, FBXSDK_CURVENODE_COMPONENT_Z);
      }

      // TODO: Weights...
    }
  }
  return true;
}

bool Raw2FbxConverter::CreateLights() {

  std::vector<FbxLight*> lights;
  for (int i = 0; i < raw.GetLightCount(); ++i) {
    auto& rawLight = raw.GetLight(i);
    FbxLight* pLight = FbxLight::Create(pScene, rawLight.name.c_str());

    switch (rawLight.type)
    {
      case RAW_LIGHT_TYPE_DIRECTIONAL:
        pLight->LightType.Set(FbxLight::eDirectional);
        break;
      case RAW_LIGHT_TYPE_POINT:
        pLight->LightType.Set(FbxLight::ePoint);
        break;
      case RAW_LIGHT_TYPE_SPOT:
        pLight->LightType.Set(FbxLight::eSpot);
        pLight->InnerAngle.Set(rawLight.innerConeAngle * 180 / M_PI);
        pLight->OuterAngle.Set(rawLight.outerConeAngle * 180 / M_PI);
        break;
    }

    pLight->Color.Set(toFbxDouble3(rawLight.color));
    pLight->Intensity.Set(rawLight.intensity);
    lights.push_back(pLight);
  }

  for (int i = 0; i < raw.GetNodeCount(); ++i) {
    // For any node that has a surface Id, add a mesh to the matching fbx node
    auto& node = raw.GetNode(i);
    if (node.lightIx != -1) {
      // Attach light to its node
      FbxNode* pLightNode = idToFbxNodeMap[node.id];
      pLightNode->SetNodeAttribute(lights[node.lightIx]);
    }
  }

  return true;
}

bool Raw2FbxConverter::CreateCameras() {
  for (int i = 0; i < raw.GetCameraCount(); ++i) {
    auto& rawCam = raw.GetCamera(i);
    FbxCamera* pCam = FbxCamera::Create(pScene, rawCam.name.c_str());

    // Set all camera properties
    switch (rawCam.mode)
    {
      case RawCamera::CAMERA_MODE_PERSPECTIVE:
        pCam->ProjectionType.Set(FbxCamera::ePerspective);
        pCam->FieldOfViewX.Set(rawCam.perspective.fovDegreesX);
        pCam->FieldOfViewY.Set(rawCam.perspective.fovDegreesY);
        pCam->FilmAspectRatio.Set(rawCam.perspective.aspectRatio);
        pCam->NearPlane.Set(rawCam.perspective.nearZ);
        pCam->FarPlane.Set(rawCam.perspective.farZ);
        break;
      case RawCamera::CAMERA_MODE_ORTHOGRAPHIC:
        pCam->ProjectionType.Set(FbxCamera::eOrthogonal);
        pCam->OrthoZoom.Set(rawCam.orthographic.magX);
        pCam->FarPlane.Set(rawCam.orthographic.farZ);
        pCam->NearPlane.Set(rawCam.orthographic.nearZ);
        break;
    }

    // Attach cam to its node
    FbxNode* pCamNode = idToFbxNodeMap[rawCam.nodeId];
    pCamNode->SetNodeAttribute(pCam);
  }

  return true;
}

bool Raw2FbxConverter::WriteNodeHierarchy() {
  
  CreateNodeHierarchy();
  
  // TODO: Write user defined properties if any
  // // Only support non-animated user defined properties for now
  // FbxProperty objectProperty = pNode->GetFirstProperty();
  // while (objectProperty.IsValid()) {
  //   if (objectProperty.GetFlag(FbxPropertyFlags::eUserDefined)) {
  //     ReadNodeProperty(raw, pNode, objectProperty);
  //   }

  //   objectProperty = pNode->GetNextProperty(objectProperty);
  // }

  CreateTextures();
  CreateMaterials();
  CreateMeshes();
  AssignMeshesAndMaterialsToNodes();
  CreateAnimations();
  CreateLights();
  CreateCameras();

  return true;
}


bool Raw2FBX(const std::string& outputPath, const RawModel& raw, const GltfOptions& options)
{
    // create a SdkManager
    FbxManager* lSdkManager = FbxManager::Create();

    // create an IOSettings object
    FbxIOSettings* ios = FbxIOSettings::Create(lSdkManager, IOSROOT);

    // set some IOSettings options 
    ios->SetBoolProp(EXP_FBX_MATERIAL, true);
    ios->SetBoolProp(EXP_FBX_TEXTURE,  true);

    // create an empty scene
    FbxScene* lScene = FbxScene::Create(lSdkManager, "");

    // glTF uses YUp, so make sure to tell the fbx scene
    lScene->GetGlobalSettings().SetAxisSystem(FbxAxisSystem::MayaYUp);
    lScene->GetGlobalSettings().SetSystemUnit(FbxSystemUnit::m);

    // FBX's internal unscaled unit is centimeters, glTF is meters, so to make things easy,
    // we pre-multiply the scale difference into every vertex position (and related attributes) instead.
    // this is always 100.0, but let's opt for clarity.
    scaleFactor = FbxSystemUnit::cm.GetConversionFactorFrom(FbxSystemUnit::m);

    // Fill it up with the raw model
    Raw2FbxConverter converter(raw, lScene);
    converter.WriteNodeHierarchy();

    // create an exporter.
    FbxExporter* lExporter = FbxExporter::Create(lSdkManager, "");

    // initialize the exporter by providing a filename and the IOSettings to use
    lExporter->Initialize(outputPath.c_str(), -1, ios);

    // export the scene.
    lExporter->Export(lScene); 

    // destroy the exporter
    lExporter->Destroy();    

    return true;
}
