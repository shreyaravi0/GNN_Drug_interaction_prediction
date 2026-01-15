# model_def.py
import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool

class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels=128):
        super(GNNEncoder, self).__init__()
        self.conv1 = GCNConv(6, hidden_channels)
        self.bn1 = torch.nn.BatchNorm1d(hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.bn2 = torch.nn.BatchNorm1d(hidden_channels)
        self.conv3 = GCNConv(hidden_channels, hidden_channels)
        self.bn3 = torch.nn.BatchNorm1d(hidden_channels)

    def forward(self, x, edge_index, edge_attr, batch):
        x = self.conv1(x, edge_index)
        x = self.bn1(x)
        x = F.relu(x)
        x = F.dropout(x, p=0.3, training=self.training)

        x = self.conv2(x, edge_index)
        x = self.bn2(x)
        x = F.relu(x)
        x = F.dropout(x, p=0.3, training=self.training)

        x = self.conv3(x, edge_index)
        x = self.bn3(x)
        x = F.relu(x)

        x = global_mean_pool(x, batch)
        return x

class InteractionPredictor(torch.nn.Module):
    def __init__(self, hidden_channels=128):
        super(InteractionPredictor, self).__init__()
        self.encoder = GNNEncoder(hidden_channels)

        # Must match the saved model naming: 'classifier.*'
        self.classifier = torch.nn.Sequential(
            torch.nn.Linear(hidden_channels * 2, hidden_channels),
            torch.nn.BatchNorm1d(hidden_channels),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.4),
            torch.nn.Linear(hidden_channels, hidden_channels // 2),
            torch.nn.BatchNorm1d(hidden_channels // 2),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.4),
            torch.nn.Linear(hidden_channels // 2, 1)
        )

    def forward(self, data1_batch, data2_batch):
        emb1 = self.encoder(data1_batch.x, data1_batch.edge_index,
                            data1_batch.edge_attr, data1_batch.batch)
        emb2 = self.encoder(data2_batch.x, data2_batch.edge_index,
                            data2_batch.edge_attr, data2_batch.batch)
        combined = torch.cat([emb1, emb2], dim=1)
        out = self.classifier(combined)
        return out
