/**
 * Copyright (c) 2019-present, Adobe, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#include "Raw2Fbx.hpp"

static double scaleFactor = 100.0; // Technically we compute it in Raw2Fbx, but it will also come out to 100 :)

static bool WriteNodeHierarchy(const RawModel& raw, FbxScene* pScene) {
  
  std::map<long, FbxNode*> idToFbxNodeMap;
  
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

  // We need a small map of surface id to node(s) that use that surface
  std::multimap<long, long> surfaceIdToNodeIdsMap;
  for (int i = 0; i < raw.GetNodeCount(); ++i) {
    // For any node that has a surface Id, add a mesh to the matching fbx node
    auto& node = raw.GetNode(i);

    // Does the node have a surface
    if (node.surfaceId >= 0) {
        // Add/Append a new entry
        surfaceIdToNodeIdsMap.insert( {node.surfaceId, node.id } );
    }
  }

  // TODO: Write user defined properties if any
  // // Only support non-animated user defined properties for now
  // FbxProperty objectProperty = pNode->GetFirstProperty();
  // while (objectProperty.IsValid()) {
  //   if (objectProperty.GetFlag(FbxPropertyFlags::eUserDefined)) {
  //     ReadNodeProperty(raw, pNode, objectProperty);
  //   }

  //   objectProperty = pNode->GetNextProperty(objectProperty);
  // }

  int attributes = raw.GetVertexAttributes();

  // Create all the meshes from surfaces
  std::map<long, FbxMesh*> surfaceIdToFbxMeshMap;

  for (int i = 0; i < raw.GetSurfaceCount(); ++i) {
    auto& surface = raw.GetSurface(i);
    FbxMesh* fbxMesh = FbxMesh::Create(pScene, surface.name.c_str());
    
    // Add the mesh to the map, so we can map it later to actual fbx nodes
    surfaceIdToFbxMeshMap[surface.id] = fbxMesh;

    // We need to collect all the triangles and verts of this surface

    // We'll store the indices of the triangles and verts used by this surface
    std::vector<int> surfaceTriIndices;
    std::vector<int> surfaceVertIndices;

    // And a mapping of the verts index in the raw model to the vert index in the collected array above
    // For this we iterate over all the triangles in the model and select the ones that match this surface.
    // With a lot of surfaces this will be very inefficient, and should instead split the process up
    // with more intermediate data and lookup.
    std::map<int, int> surfaceVertIndexToMeshVertIndex;
    for (int t = 0; t < raw.GetTriangleCount(); ++t)
    {
      auto& tri = raw.GetTriangle(t);
      if (tri.surfaceIndex == i) {
        surfaceTriIndices.push_back(t);
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

    // RAW_VERTEX_ATTRIBUTE_JOINT_INDICES = 1 << 7,
    // RAW_VERTEX_ATTRIBUTE_JOINT_WEIGHTS = 1 << 8,

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

    // Skinning
    if (attributes & RAW_VERTEX_ATTRIBUTE_JOINT_INDICES || attributes & RAW_VERTEX_ATTRIBUTE_JOINT_WEIGHTS) {

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

      // Create the skeletons
      auto& rawSkeletonRootNode = raw.GetNode(raw.GetNodeById(surface.skeletonRootId));
      if (!rawSkeletonRootNode.isJoint) {
        continue;
      }

      FbxSkeleton* pSkeletonRootAttribute = FbxSkeleton::Create(pScene, rawSkeletonRootNode.name.c_str());
      pSkeletonRootAttribute->SetSkeletonType(FbxSkeleton::eLimb);
      FbxNode* pSkeletonRoot = idToFbxNodeMap[surface.skeletonRootId];
      pSkeletonRoot->SetNodeAttribute(pSkeletonRootAttribute);

      // Now iterate through the children and children of children
      //

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
      fbxMesh->AddDeformer(pSkin);
    }
  }

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
    }
  }

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
    WriteNodeHierarchy(raw, lScene);

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
